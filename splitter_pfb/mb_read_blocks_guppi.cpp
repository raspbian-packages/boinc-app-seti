#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <vector>
#include "fitsio.h"
#include "psrfits.h"
#include "guppi_params.h"
#include "fitshead.h"
#include "median.h"
#include "setimysql.h"
#include <fftw3.h>
#include <sys/stat.h>
#include "barycenter.h"
#include "rawdopplersearch.h"
#include "setilib.h"
#include "mb_read_blocks_guppi.h"
#include "mb_splitter.h"
#include "mb_coords.h"

#define FIX_PARKES_HEADER

//#include "seti_dr2utils.h"

/* Guppi channel-frequency mapping */
/* sample at 1600 MSamp for an 800 MHz BW */
/* N = 256 channels with frequency centers at m * fs/N */
/* m = 0,1,2,3,... 255 */

// guppi raw data layout (data portion of a ~psrfits HDU (subintegration)):
//
// First of 32 frequency channels / subintegration:
// Pol0Samp0    Pol0Samp0   Pol1Samp0   Pol1Samp0   ... Pol1Samp255 Pol1Samp255
// real         imag        real        imag            real        imag
// |--------------------------------------------|       ----------------------|
//   byte 0 (2 bit data, quantized from 8 bit)              ... byte 255
//   time 0                                                 ... time 255
//
// The next subintegration will be the same 32 channels, starting with time 256.
// Thus, each subintegration is 32*256 = 8192 bytes of data.
//
#if 0

We read the actual data samples ass pointed to by rawinput.pf.sub.data

struct gpu_input {                                                              // this is instantiated as rawdata
    char *file_prefix;
    struct                  guppi_params gf;
    struct                  psrfits pf;                                         // path to data : rawdata.pf
    unsigned int filecnt;
    FILE *fil;
    int invalid;
    int curfile;
    int overlap;   /* add this keyword here since it doesn't seem to appear in guppi_params.c */
    long int first_file_skip; /* in case there's 8bit data in the header of file 0 */
};

struct guppi_params {
    /* Packet information for the current block */
    long long packetindex;      // Index of first packet in raw data block
    int packetsize;             // Size in bytes of data portion of each packet
    int n_packets;              // Total number of packets in current block
    int n_dropped;              // Number of packets dropped in current block
};

struct psrfits {
    char mode;              // Read (r) or write (w).
    struct                  hdrinfo hdr;
    struct                  subint sub;                                         // path to data : rawdata.pf.sub
};

struct hdrinfo {
    char source[24];        // Source name
    char frontend[24];      // Frontend used
    long double MJD_epoch;  // Starting epoch in MJD
    double fctr;            // Center frequency of the observing band (MHz)
    double df;              // Frequency spacing between the channels (MHz)		// in BL data, this is header item CHAN_BW
    double BW;              // Bandwidth of the observing band (MHz)			// in BL data, this is header item OBSBW
    int nbits;              // Number of bits per data sample
    int nchan;              // Number of channels
    int rcvr_polns;         // Number of polns provided by the receiver
};

struct subint {
    double tsubint;         // Length of subintegration (sec)
    double offs;            // Offset from Start of subint centre (sec)
    double ra;              // RA (J2000) at subint centre (deg)
    double dec;             // Dec (J2000) at subint centre (deg)
    int bytes_per_subint;   // Number of bytes for one row of raw data
    unsigned char *data;    // Ptr to the raw data itself                       // path to data : rawdata.pf.sub.data
};

#endif

unsigned long total_samples = 0;
unsigned long limit_samples;
unsigned long last_num_sample_bytes_read;

#define RAW_DATA_HEADER_BUF_SIZE 32768
//#define RAW_DATA_HEADER_BUF_SIZE 8192         // TODO would it help to do less reading/rewinding for each header?

typedef struct {
    int             guppi_type;
    // file control...
    long long int   startindx;
    long long int   curindx;
	long long int   chanbytes;      // = control_block.indxstep * rawinput.gf.packetsize / (8/raw_bit_depth * rawinput.pf.hdr.nchan); 
	long long int   subint_offset;
	long int        currentblock;
	int             indxstep;       // = (int) ((rawinput.pf.sub.bytes_per_subint * samples_per_byte) / rawinput.gf.packetsize) - 
                                    //   (int) (rawinput.overlap * rawinput.pf.hdr.nchan * rawinput.pf.hdr.rcvr_polns           * 
                                    //   2 / rawinput.gf.packetsize);
    int             last_block_done;
    // user prefs...
    int             channel;
    int             polarization;
    long int        startblock; 
    long int        numblocks; 
    int             vflag; 
} control_block_t;

/* Parse info from buffer into param struct */
extern void guppi_read_obs_params(char *buf, 
                                     struct guppi_params *g,
                                     struct psrfits *p);

/* prototypes */

int exists(const char *fname);

unsigned long unpack_samples(unsigned char * raw, long int count, int pol, dr2_compact_block_t &this_tapebuffer);

//-------------------------------------------------------
void print_header(dr2_compact_block_t &this_tapebuffer) {
//-------------------------------------------------------

	fprintf(stderr, "setiheader : Name: %s Frontend: %s Header Size: %d Data Size: %d Channel: %d Polarization: N samplerate %lf sky_freq %ld Data Time: %lf Coordinate Time: %lf RA: %f Dec: %f\n", 
            this_tapebuffer.header.name,
            this_tapebuffer.header.frontend,
            this_tapebuffer.header.header_size,
	        this_tapebuffer.header.data_size,		
	        this_tapebuffer.header.channel,		
	        this_tapebuffer.header.samplerate,		
	        this_tapebuffer.header.sky_freq,		
            this_tapebuffer.header.data_time.mjd().uval(),
            this_tapebuffer.header.coord_time.mjd().uval(),
	        this_tapebuffer.header.ra,
	        this_tapebuffer.header.dec);
}

//-------------------------------------------------------
int get_frontend(struct gpu_input &rawinput, dr2_compact_block_t &this_tapebuffer) {
//-------------------------------------------------------

	int rv = 1;		// assume no error

    strncpy(this_tapebuffer.header.frontend, rawinput.pf.hdr.frontend, sizeof(this_tapebuffer.header.frontend));
	
	if(!strncmp(this_tapebuffer.header.frontend, "1050CM", strlen("1050CM"))) {
		// the Parkes 1050CM is dual band. We differentiate here.
		if(rawinput.pf.hdr.fctr >= 2600.0 && rawinput.pf.hdr.fctr <= 3600.0) {			// between 2.6 and 3.6 GHz
			strcat(this_tapebuffer.header.frontend, "_10"); 
		} else if (rawinput.pf.hdr.fctr >= 700.0 &&  rawinput.pf.hdr.fctr <= 764.0) {	// between 700 and 764 MHz
			strcat(this_tapebuffer.header.frontend, "_50"); 
		} else {
			fprintf(stderr, "Frequency %lf is inconsistent with frontend %s\n",  
					rawinput.pf.hdr.fctr, rawinput.pf.hdr.frontend);
			rv = 0;
		}
	} else if(!strncmp(this_tapebuffer.header.frontend, "BL", strlen("BL"))) {
		char beam_id[16];
		sprintf(beam_id, "%d", rawinput.pf.hdr.beam_id);		
		// this is Parkes multibeam data and so we must add the beam ID
		// to an otherwise ambiguous frontend name.  In fact, we change
		// the frontend name completely for the sake of clarity.
		strcpy(this_tapebuffer.header.frontend, "Parkes_MB_");
		strcat(this_tapebuffer.header.frontend, beam_id);
	}

			
	return(rv);
}

//-------------------------------------------------------
int populate_header(
                char                *file_prefix,
	            control_block_t     &control_block,
                struct gpu_input    &rawinput,
                dr2_compact_block_t &this_tapebuffer
) {
//-------------------------------------------------------

	int rv = 1; 	// assume no error

#if 0
	// orig frequency code:
    // this assumes that channels are ordered with *descending* frequency 
    // (the higher the channel, the lower the frequency) and that BW
    // and df are always positive.
    // TODO - channel ordering should come from a recorder_config field
    this_tapebuffer.header.sky_freq =   (long int) (
                              (double) 1e6 * 
                              ( 
                                (rawinput.pf.hdr.fctr   + 
                                (rawinput.pf.hdr.BW/2)) -      
                                ((control_block.channel+0.5) * rawinput.pf.hdr.df)
                              )
                            );
#endif
//#if 0
	// new frequency code:
    // this assumes in the case of *descending* frequency, both BW and df will
    // be negative, otherwise positive.  Ie, this should work in all cases where
    // the -BW and -df protocol is followed. 
	this_tapebuffer.header.sky_freq =   (long int) ( 
							(double) 1e6 *  
	                        (  
 	                          (rawinput.pf.hdr.fctr   -  
	 	                      (rawinput.pf.hdr.BW/2)) +  
	                          ((control_block.channel+0.5) * rawinput.pf.hdr.df) 
	                        ) 
						  );
//#endif

fprintf(stderr, "rawinput.pf.hdr.BW is %f\n", rawinput.pf.hdr.BW);
fprintf(stderr, "rawinput.pf.hdr.df is %f\n", rawinput.pf.hdr.df);
	if(rawinput.pf.hdr.BW < 0) {
 		log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG,"Adjusting for negative BW : rawinput.pf.hdr.BW was %f ... ", rawinput.pf.hdr.BW); 
		rawinput.pf.hdr.BW = fabs(rawinput.pf.hdr.BW);
 		log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG,"is %f ", rawinput.pf.hdr.BW); 
	}
	if(rawinput.pf.hdr.df < 0) {
		log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG, "Adjusting for negative BW : rawinput.pf.hdr.df was %f ... ", rawinput.pf.hdr.df); 
		rawinput.pf.hdr.df = fabs(rawinput.pf.hdr.df);
		log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG,"is %f ... ", rawinput.pf.hdr.df); 
	}
fprintf(stderr, "rawinput.pf.hdr.BW is %f\n", rawinput.pf.hdr.BW);
fprintf(stderr, "rawinput.pf.hdr.df is %f\n", rawinput.pf.hdr.df);

    /* size of this header */
    this_tapebuffer.header.header_size = sizeof(this_tapebuffer.header);

    /* populate with source name from psrfits header */
    strncpy(this_tapebuffer.header.name, basename(file_prefix), sizeof(this_tapebuffer.header.name));

    /* populate with frontend name from psrfits header */
    //strncpy(this_tapebuffer.header.frontend, rawinput.pf.hdr.frontend, sizeof(this_tapebuffer.header.frontend));
	rv = get_frontend(rawinput, this_tapebuffer);

    /* populate channel and polarization */
    this_tapebuffer.header.channel      = control_block.channel;
    // this_tapebuffer.header.polarization = control_block.polarization;   // TODO add polarization

    /* (complex) sample rate in MHz */
    this_tapebuffer.header.samplerate = rawinput.pf.hdr.df;

    /* there is no "dataseq in guppi data.  TODO - right? */
    //this_tapebuffer.header.dataseq = 0;
    print_header(this_tapebuffer);

	return(rv);
}


//-------------------------------------------------------
int get_params(
	            control_block_t     &control_block,
                struct gpu_input    &rawinput,
                char                header_buf[],
                char                filname[]
) {
//-------------------------------------------------------
	long long int chanbytes_overlap = 0;
        // Read first raw data header
	    if(fread(header_buf, sizeof(char), RAW_DATA_HEADER_BUF_SIZE, rawinput.fil) == RAW_DATA_HEADER_BUF_SIZE){
		   
		    guppi_read_obs_params(header_buf, &rawinput.gf, &rawinput.pf);
#ifdef FIX_PARKES_HEADER
			if(control_block.guppi_type == guppi_parkes) {
				rawinput.pf.hdr.nbits = 2;
				rawinput.pf.sub.bytes_per_subint /= 4;
			}
#endif
            // Get to first N bit data header where N != 8...
            // At regular intervals, as a sanity check, the N bit
            // data are preceeded by the full 8 bit version of
            // the same data.  We ignore the 8 bit data here.
            // We ignore this whole check if we have 8 bit data.
		    if(raw_bit_depth != 8 && rawinput.pf.hdr.nbits == 8) {
			    fprintf(stderr, "caught an 8 bit header... moving on to first/next 2 bit header\n");
			    /* figure out the size of the first subint + header */
		        rawinput.first_file_skip = rawinput.pf.sub.bytes_per_subint + gethlength(header_buf);
			    /* rewind to the beginning */	
		   	    fseek(rawinput.fil, -RAW_DATA_HEADER_BUF_SIZE, SEEK_CUR);
		   	    /* seek past the first subint + header */
		   	    fseek(rawinput.fil, rawinput.first_file_skip, SEEK_CUR);
			    /* read the next header */
		   	    fread(header_buf, sizeof(char), RAW_DATA_HEADER_BUF_SIZE, rawinput.fil);
			    guppi_read_obs_params(header_buf, &rawinput.gf, &rawinput.pf);
			    fclose(rawinput.fil);
		    } else {
                // not 8 bit data 
			    fclose(rawinput.fil);
		    }   // end if(rawinput.pf.hdr.nbits == 8)

		    /* we'll use this header to set the params for the whole observation */
		   
		    rawinput.fil = NULL;

		    hgeti4(header_buf, "OVERLAP", &rawinput.overlap);
			
		    fprintf(stderr, "   pktindx           : %lld\n", rawinput.gf.packetindex);
		    fprintf(stderr, "   packetsize        : %d\n", rawinput.gf.packetsize);
		    fprintf(stderr, "   n_packets         : %d\n", rawinput.gf.n_packets);
		    fprintf(stderr, "   n_dropped         : %d\n",rawinput.gf.n_dropped);
		    fprintf(stderr, "   bytes_per_subint  : %d\n",rawinput.pf.sub.bytes_per_subint);
		    fprintf(stderr, "   overlap           : %d\n",rawinput.overlap);
		    fprintf(stderr, "   bit depth         : %d\n",raw_bit_depth);
		    fprintf(stderr, "   beam_id           : %d\n", rawinput.pf.hdr.beam_id);
	        } else {    
	  		    fprintf(stderr, "couldn't read a header\n");
			    return 0;
	        }   // end fread(buf...


    if(control_block.vflag>=1) fprintf(stderr, "calculating index step\n");

    /* number of packets that we *should* increment by */
    //control_block.indxstep = (int) ((rawinput.pf.sub.bytes_per_subint * 4) / rawinput.gf.packetsize) - 
    control_block.indxstep = (int) ((rawinput.pf.sub.bytes_per_subint * samples_per_byte) / rawinput.gf.packetsize) - 
                      (int) (rawinput.overlap * rawinput.pf.hdr.nchan * rawinput.pf.hdr.rcvr_polns  * 
                      2 / rawinput.gf.packetsize);

    if(control_block.vflag>=1) {
        fprintf(stderr, "index step %d bytes_per_subint %d packetsize %d overlap %d nchan %d rcr_polns %d\n",
                control_block.indxstep, rawinput.pf.sub.bytes_per_subint, rawinput.gf.packetsize,
                rawinput.overlap, rawinput.pf.hdr.nchan, rawinput.pf.hdr.rcvr_polns);
    }

    if (control_block.channel >= rawinput.pf.hdr.nchan) {
	    fprintf(stderr, "control_block.channel %d more than channels in data %d\n", control_block.channel, rawinput.pf.hdr.nchan);
	    return 0;
    } else {
	    fprintf(stderr, "Numer of channels in file %d\n", rawinput.pf.hdr.nchan);
    }

    /* number of non-overlapping bytes in each channel */
    /* control_block.indxstep increments by the number of unique packets in each sub-integration */
    /* packetsize is computed based on the original 8 bit resolution */
    /* divide by 8/raw_bit_depth to get to 2 or 8 bits, nchan to get to number of channels */
    // TODO - this needs to change for 8 bit data!!
    //control_block.chanbytes = control_block.indxstep * rawinput.gf.packetsize / (4 * rawinput.pf.hdr.nchan); 
    control_block.chanbytes = control_block.indxstep * rawinput.gf.packetsize / (8/raw_bit_depth * rawinput.pf.hdr.nchan); 
    fprintf(stderr, "chan bytes %lld with a bit depth of %d\n", control_block.chanbytes, raw_bit_depth);

    /* total number of bytes per channel, including overlap */
    chanbytes_overlap = rawinput.pf.sub.bytes_per_subint / rawinput.pf.hdr.nchan;

    /* memory offset for our chosen channel within a subint */
    control_block.subint_offset = control_block.channel * chanbytes_overlap;

    if(control_block.vflag>=1) fprintf(stderr, "Index step: %d\n", control_block.indxstep);
    if(control_block.vflag>=1) fprintf(stderr, "bytes per subint %d\n",rawinput.pf.sub.bytes_per_subint );
}


//-------------------------------------------------------
int read_block(
	            control_block_t     &control_block,
                struct gpu_input    &rawinput,
                char                header_buf[],
                char                filname[],
                double              &start_data_time,
                double              &end_data_time,
                long int            filecnt,
                dr2_compact_block_t &this_tapebuffer
) {
//-------------------------------------------------------

	size_t rv=0;
    int retval=0;
	//long int j=0;

	    if(!rawinput.invalid){						  
		    if(rawinput.fil == NULL) {
			    /* no file is open for this band, try to open one */
			    sprintf(filname, "%s.%04d.raw",rawinput.file_prefix,rawinput.curfile);
			    if(exists(filname)){
                    fprintf(stderr, "Opening %s\n", filname);
				    rawinput.fil = fopen(filname, "rb");			 
				    if(rawinput.curfile == 0 && rawinput.first_file_skip != 0) fseek(rawinput.fil, rawinput.first_file_skip, SEEK_CUR);  
			    }	else {
			  	    rawinput.invalid = 1;
		  	  	    fprintf(stderr, "couldn't open any more files!\n");
		  	    }
		    }  // end  if(rawinput.fil == NULL)

		    if(rawinput.fil){		  
		        fprintf(stderr, "Reading %ld bytes - more than enough to obtain first/next header (we'll rewind)\n", RAW_DATA_HEADER_BUF_SIZE); 
	            if(fread(header_buf, sizeof(char), RAW_DATA_HEADER_BUF_SIZE, rawinput.fil) == RAW_DATA_HEADER_BUF_SIZE) {
		            //fprintf(stderr, "...success!\n"); 
					fseek(rawinput.fil, -RAW_DATA_HEADER_BUF_SIZE, SEEK_CUR);           // rewind for (header plus) data read

					if(control_block.vflag>=1) fprintf(stderr, "header length: %d\n", gethlength(header_buf));

                    // parse the header we just read
					guppi_read_obs_params(header_buf, &rawinput.gf, &rawinput.pf);
#ifdef FIX_PARKES_HEADER
					if(control_block.guppi_type == guppi_parkes) {
						rawinput.pf.hdr.nbits = 2;
						rawinput.pf.sub.bytes_per_subint /= 4;
					}
#endif

                    // get our current block number 
					control_block.currentblock = (long int) ((double) rawinput.gf.packetindex/ (double) control_block.indxstep);

					if(control_block.vflag>=1) {
						 fprintf(stderr, "pktindx %Ld packetsize: %d n_packets %d n_dropped: %d blocks: %ld RA: %f DEC: %f subintoffset %f tsubint %f MJD %Lf\n", 
                                    rawinput.gf.packetindex,
						            rawinput.gf.packetsize,
						            rawinput.gf.n_packets,
						            rawinput.gf.n_dropped,
						            control_block.currentblock,
						            rawinput.pf.sub.ra,
						            rawinput.pf.sub.dec,
						            rawinput.pf.sub.offs,
						            rawinput.pf.sub.tsubint,
						            rawinput.pf.hdr.MJD_epoch);
					}

					
					/* populate the variables that change with each block */
				    /* RA (J2000) at subint centre (deg), Dec (J2000) at subint centre (deg), time in MJD */
                    this_tapebuffer.header.dataseq      = control_block.currentblock;
					this_tapebuffer.header.ra           = rawinput.pf.sub.ra/15.0;      // splitter expects RA in decimal hours            
					this_tapebuffer.header.dec          = rawinput.pf.sub.dec;           
					/* increment time for data by 0.5 x length of block to push quoted time to the _last_ sample in the block */					
					double data_mjd                           = rawinput.pf.hdr.MJD_epoch + ((rawinput.pf.sub.offs + rawinput.pf.sub.tsubint/2)/86400.0);     
					double coord_mjd                          = rawinput.pf.hdr.MJD_epoch + (rawinput.pf.sub.offs/86400.0);    
                    this_tapebuffer.header.data_time    = seti_time(days(data_mjd),MJD0);
                    this_tapebuffer.header.coord_time   = seti_time(days(coord_mjd),MJD0);

                    if(control_block.currentblock == 1) {
                        start_data_time = this_tapebuffer.header.data_time.mjd().uval();
                    } else {
                        end_data_time = this_tapebuffer.header.data_time.mjd().uval();     // assigned many times, will end up with the final block's time
                    }
					
			   		if(rawinput.gf.packetindex == control_block.curindx) {    // file integrity check - are we at the correct index?  Yes...

						 long seek_len;
						 long hlength = (long)gethlength(header_buf);
						 if(rawinput.pf.hdr.directio) {
							seek_len = hlength + (512 - (hlength%512))%512;				// directio aligns on 512 byte boundaries
							fprintf(stderr, "adjusting for directio (%d) : %ld becomes %ld\n", rawinput.pf.hdr.directio, hlength, seek_len);
						 } else {
							seek_len = hlength;											// not directio so no forced alignment
							fprintf(stderr, "NOT adjusting for directio (%d) : %ld becomes %ld\n", rawinput.pf.hdr.directio, hlength, seek_len);
						 }
						 fseek(rawinput.fil, seek_len, SEEK_CUR);         				  // skip past header
						 rv=0;
						 
						 if (control_block.currentblock > (control_block.startblock-1)){   // if we have reached user requested startblock (default 0)
						 
                            // read raw data.  bytes_per_subint should == BLOCSIZE value from header
						 	rv=fread(rawinput.pf.sub.data, sizeof(char), rawinput.pf.sub.bytes_per_subint, rawinput.fil); 
						 	
							if( ((long int)rv == rawinput.pf.sub.bytes_per_subint) ){       // if we got the expected amount of raw data
                               long filepos = ftell(rawinput.fil);							   
                               if(control_block.vflag>=1) {
                                    fprintf(stderr,"read %d bytes starting at filepos %ld in file %ld of %d ... chanbytes is %d\n", 
                                            rawinput.pf.sub.bytes_per_subint, filepos, rawinput.curfile, filecnt, control_block.chanbytes);						
                                    }
                               // wrap unpacking samples in restart logic.  We skip the unpack (but do everything else) if we are not
                               // to the point of resumption.
                               if(control_block.currentblock >= control_block.last_block_done) {
                                    //if(control_block.vflag>=1) fprintf(stderr, "prior to unpack_samples : last : %ld total : %ld\n", 
                                    //                                   last_num_sample_bytes_read, total_samples);
						            last_num_sample_bytes_read = unpack_samples(rawinput.pf.sub.data     +   // start of raw data, ie channel 0
                                                                           control_block.subint_offset, // this gets us to the requested channel 
                                                                           control_block.chanbytes,     // chanbytes does not include overlap
                                                                           control_block.polarization,
                                                                           this_tapebuffer);   			// unpack samples to the return vector   
                                    if(control_block.vflag>=1) fprintf(stderr, "Just unpacked %ld sample bytes from %ld chanbytes.  Total samples : %ld\n", 
                                                                       last_num_sample_bytes_read, control_block.chanbytes, total_samples);
                               } else {
                                    log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG,"Not to point of resumption (%ld < %ld), skipping...\n", 
                                                        control_block.currentblock, control_block.last_block_done);
                               }  // end of unpack with restart logic

							} else {
								rawinput.fil = NULL;
								rawinput.invalid = 1;
								fprintf(stderr,"ERR: couldn't read as much as the header said we could... assuming corruption and exiting...\n");
								exit(1);
							}  // end if we got the expected amount of raw data

						 } else {                                                                   // we have not reached the requested startblock
                            (rv = fseek(rawinput.fil, rawinput.pf.sub.bytes_per_subint, SEEK_CUR)); // ... so seek to the next block
                         }  // end if we have reached user requested startblock
						 
					} else if( (rawinput.gf.packetindex > control_block.curindx) && (control_block.currentblock > (control_block.startblock-1)) ) {
                         // are we at the correct index?  No, we apparently dropped a subintegration.  This is correctable...

						 fprintf(stderr,"ERR: control_block.curindx: %Ld, pktindx: %Ld Did we drop a whole subintegration? Reusing last subint...\n", 
                                 control_block.curindx, rawinput.gf.packetindex );
						 /* read a subint with too high an indx, must have dropped a whole subintegration*/

						/* pf.sub.data *should* still contain the last valid subint */
						/* grab a copy of the last subint  - probably should add gaussian noise here, but this is better than nothing! */

						/* we'll keep the ra/dec values from this subint, but push back the time to keep everything sensible */
					    data_mjd = rawinput.pf.hdr.MJD_epoch + ((rawinput.pf.sub.offs - rawinput.pf.sub.tsubint/2)/86400.0) ;     
                        this_tapebuffer.header.data_time = seti_time(days(data_mjd),MJD0);

                        // TODO keep the following demo buffer code to remind me what to do here
                        // header to return buffer
						//memcpy((*setibuffer) + control_block.setibuffer_pos, 
                        //       &setiheader, 
                        //       sizeof(setiheader));
						//control_block.setibuffer_pos += sizeof(setiheader);

                        // data to return buffer
						//memmove((*setibuffer) + control_block.setibuffer_pos, 
                        //        (*setibuffer) + control_block.setibuffer_pos - sizeof(setiheader) - control_block.chanbytes, 
                        //        (control_block.chanbytes * 2 * sizeof(unsigned char)));
						//control_block.setibuffer_pos += (control_block.chanbytes * 2 * sizeof(unsigned char));

						/* We'll get the current valid subintegration again on the next time through this loop */

					} else if(rawinput.gf.packetindex < control_block.curindx) {
                         // are we at the correct index?  No, nor have we dropped a subintegration.  This is a file integrity fatal error...
						 fprintf(stderr,"Error expecting a higher packet index than we got control_block.curindx: %Ld, pktindx: %Ld\n", 
                                 control_block.curindx, rawinput.gf.packetindex );
						 /* somehow we were expecting a higher packet index than we got !?!? */	
						 fprintf(stderr, "assuming corruption and exiting...\n");
						 exit(1);
					}  // end are we at the correct index
    
                if(control_block.vflag>=1) fprintf(stderr, "Read block %ld\n", control_block.currentblock);

                //if(control_block.vflag>=1) print_header(this_tapebuffer);
                retval = true;
				} else {
				   fprintf(stderr, "Could not read first/next 32KB\n"); 
				   fclose(rawinput.fil);
				   rawinput.fil = NULL;
				   rawinput.curfile++;						
                   retval = false;
				}   // end if((fread(buf...)
		    }   // end if(rawinput.fil)			 	 	 
	    }   // end if(!rawinput.invalid)
										
	    if(rawinput.fil != NULL) {
            //fprintf(stderr, "index step : %d\n", control_block.indxstep);
            control_block.curindx = control_block.curindx + control_block.indxstep;
            retval = 1;     // we have a new block, which may be empty if we are fast forwarding to the point of resumption
        } else {
            retval = 0;     // we do not have a new block. Perhaps we just need to open the next file in the series.
        }

    //if(control_block.vflag>=1) fprintf(stderr, "last_num_sample_bytes_read %ld chanbytes %ld\n", last_num_sample_bytes_read, control_block.chanbytes);
    return(retval);
}

//-------------------------------------------------------
int init(
                int                 guppi_type,
                char *              file_prefix,
                long int            startblock, 
                int                 channel, 
                int                 polarization, 
                int                 vflag, 
	            control_block_t     &control_block,
                struct gpu_input    &rawinput,
                char                header_buf[],
                char                filname[],
                double              &start_data_time,
                double              &end_data_time,
                int                 &filecnt,
                dr2_compact_block_t &this_tapebuffer
) {
//-------------------------------------------------------

    rawinput.file_prefix        =file_prefix;
    rawinput.fil                = NULL;
    rawinput.invalid            = 0;
	rawinput.first_file_skip    = 0;  
    //control_block.check_last    = 1;
    control_block.guppi_type    = guppi_type;
    control_block.chanbytes     = 0;
    control_block.subint_offset = 0;
    control_block.currentblock  = 0;
    control_block.indxstep      = 0;
    control_block.channel       = channel;
    control_block.polarization  = polarization;
    //control_block.numblocks     = numblocks;
    control_block.startblock    = startblock;
    control_block.vflag         = vflag;

    if(rawinput.file_prefix == NULL) {
	    printf("ERR no input stem specified, exiting...\n");
	    exit(1);
    }

    if(strstr(rawinput.file_prefix, ".0000.raw") != NULL) memset(rawinput.file_prefix + strlen(rawinput.file_prefix) - 9, 0x0, 9);

    // get size of data set, in terms of file count and byte count
    fprintf(stderr, "Finding size of data set...\n");
    //j = 0;
    struct stat st;
    long int size=0;
	long int i=0;
    do {
	    sprintf(filname, "%s.%04ld.raw",rawinput.file_prefix,i);
	    fprintf(stderr, "  trying %s ...",filname);		
	    i++;
	    if(exists(filname)) { 
            fprintf(stderr, " found\n");
		    stat(filname, &st);
		    size = size + st.st_size;
	    } else {
            fprintf(stderr, " not found\n");
        }
    } while (exists(filname));
    rawinput.filecnt = i-1;
    fprintf(stderr, "  File count is %i  Total data set size is %ld bytes\n",rawinput.filecnt, size);

   /* didn't find any files */
    if(rawinput.filecnt < 1) {
	    fprintf(stderr, "no files for stem %s found, exiting...\n",rawinput.file_prefix);
	    exit(1);		
    }

    /* open the first file for input */
    sprintf(filname, "%s.0000.raw", rawinput.file_prefix);
    fprintf(stderr, "Opening %s\n", filname);
    rawinput.fil = fopen(filname, "rb");
    if(!rawinput.fil){
        fprintf(stderr, "couldn't open first file\n");
	    return 0;
    }  

    get_params(
                control_block,
                rawinput,
                header_buf,
                filname
    );

    // Now that we know our sizing requirements, allocate the buffers...
    // TODO - we should explicitly free this memory even though we allocate
    // once and need it until program exit.
    rawinput.pf.sub.data  = (unsigned char *) malloc(rawinput.pf.sub.bytes_per_subint);
    if(!rawinput.pf.sub.data) {
	    fprintf(stderr, "error: couldn't allocate memory for the raw data read buffer\n");
	    return 0;
    } else {
        fprintf(stderr, "malloc'ed %ld bytes at %p for the raw data read buffer!\n", rawinput.pf.sub.bytes_per_subint, rawinput.pf.sub.data);
    }

    // populate the header with static items
    if(!populate_header(
                rawinput.file_prefix,
	            control_block,
                rawinput,
                this_tapebuffer
    )) {
		return 0;	// return error
	}
    
    control_block.startindx        = rawinput.gf.packetindex;
    control_block.curindx          = control_block.startindx;
    filecnt                        = rawinput.filecnt;
    rawinput.curfile               = 0;			

    return 1;   // success
}

//-------------------------------------------------------
void dump_tapebuffer(char * file_prefix, int channel) {
//-------------------------------------------------------

        FILE * dump_fp;
        char dump_fn[256];
        sprintf(dump_fn, "%s%s%d", basename(file_prefix), ".tapebuffer_dump.channel_", channel);
        dump_fp = fopen(dump_fn, "w");
        if(dump_fp == NULL) {
             fprintf(stderr, "could not open dump file %s, exiting...", dump_fn);
                exit(1);
        }
        for(int buffer_i=0; buffer_i < tapebuffer.size(); buffer_i++) {
//fprintf(stderr, "buffer_i %d\n", buffer_i);
            for(int buffer_j=0; buffer_j < tapebuffer[buffer_i].data.size(); buffer_j++) {
//fprintf(stderr, "buffer_i %d buffer_j %ld\n", buffer_i, buffer_j);
                fwrite((const void *)&tapebuffer[buffer_i].data[buffer_j].real(), 1, 1, dump_fp);
                fwrite((const void *)&tapebuffer[buffer_i].data[buffer_j].imag(), 1, 1, dump_fp);
            }
        }
        fclose(dump_fp);
}



//-------------------------------------------------------
long int read_blocks_guppi(char *           file_prefix,
                           int              guppi_type,
                           long int         startblock, 
                           long int         numsamples, 
                           int              channel, 
                           int              polarization, 
                           int              vflag 
) {
//-------------------------------------------------------
// Note : a buffer is push_back()'ed to tapebuffer under two conditions:
// 1) at the end of processing a guppi block.  This is by far the most
//      common push_back() and will result in a tapebuffer of a standard
//      length for a given file format until we reach the desired number
//      of samples per call. The final buffer will likely be a partial
//      buffer.
// 2) Upon re-entry into read_blocks_guppi() for consecutive WUGs, we 
//      push_back() a buffer with the left-over samples in the prior guppi 
//      block. This will be a partial buffer.  The number of samples in 
//      this buffer plus the number of samples in the final (partial) buffer 
//      of the prior call should be the standard length.      

    int retval = 0;
    static bool first_time = true;

    char header_buf[RAW_DATA_HEADER_BUF_SIZE];   // we read the raw data headers into this buffer

    long reent_total_samples = 0;

    static unsigned long sample_bytes_left;
    static double start_data_time;
    static double end_data_time;
	static int filecnt;
	static char filname[250];
    static dr2_compact_block_t this_tapebuffer;
	static struct gpu_input rawinput;	
    static control_block_t control_block;

    limit_samples = numsamples;                 // TODO - this should go into the control block
    total_samples = 0;
     
    if(first_time) {
        first_time = false;
        retval = init(
                    guppi_type,
                    file_prefix,
                    startblock,
                    channel,
                    polarization,
                    vflag,                
                    control_block,
                    rawinput,
                    header_buf,
                    filname,
                    start_data_time,
                    end_data_time,
                    filecnt,
                    this_tapebuffer
                 );
        if(retval == 0) return(0);
        if (resumetape) {
            char * filename = basename(rawinput.file_prefix);
            tape thistape;
            thistape.id=0; 
            char beambuf[16];
            sprintf(beambuf,"%d", 2*channel+polarization);            
            if (thistape.fetch(std::string("where name=\'")+filename+"\' and beam="+beambuf)) {
                control_block.last_block_done = thistape.last_block_done;
                log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG,
                                    "Resuming tape %s channel %d pol %d at block %d\n",
                                    thistape.name,channel,polarization,control_block.last_block_done);
            } else {
                log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG,
                                    "Cannot resume : tape %s channel %d pol %d beam %s not found.  Starting at the beginning.\n",
                                    filename,channel,polarization, beambuf);
                //return(0);
            }
        }  // end if(resumetape)
    } else {    // not first time -  on subsequent calls we need to account for splitter (not guppi) overlap
        reent_total_samples = 0;
        for(int buffer_i=0; buffer_i < tapebuffer.size(); buffer_i++) {
            reent_total_samples += tapebuffer[buffer_i].data.size();
            if(control_block.vflag>=2) {
                fprintf(stderr, "re-entry (start) - tapebuffer %d has %ld samples\n", buffer_i, tapebuffer[buffer_i].data.size());
            }
        }
        total_samples = reent_total_samples;        // starting this call with thais many samples :
        if(control_block.vflag>=2) {
            fprintf(stderr, "re-entry (start) - total_samples = %ld (%ld) in %ld tapebuffers\n", reent_total_samples, total_samples, tapebuffer.size());
        }
    }  // end first/subsequent time logic 


    // take care of any samples left over in the guppi buffer
    if(sample_bytes_left) {
        if(control_block.vflag>=1) fprintf(stderr, "first getting the %ld sample_bytes_left starting at location %p\n", 
                                           sample_bytes_left, rawinput.pf.sub.data + control_block.subint_offset + control_block.chanbytes - sample_bytes_left);
	    last_num_sample_bytes_read = unpack_samples(rawinput.pf.sub.data + control_block.subint_offset + control_block.chanbytes - sample_bytes_left, 
                                              sample_bytes_left,                          // samples left does not include overlap
                                              control_block.polarization,
                                              this_tapebuffer);   				    // unpack samples to the return vector   
        tapebuffer.push_back(this_tapebuffer);                                      // this will be a partially full buffer
        this_tapebuffer.data.clear();       
        sample_bytes_left = 0;
        if(control_block.vflag>=2) {
            long reent_total_samples = 0;
            for(int buffer_i=0; buffer_i < tapebuffer.size(); buffer_i++) {
                reent_total_samples += tapebuffer[buffer_i].data.size();
                fprintf(stderr, "re-entry (after adding leftover) tapebuffer %d has %ld samples\n", buffer_i, tapebuffer[buffer_i].data.size());
            }
            fprintf(stderr, "re-entry (after adding leftover) - total_samples = %ld (%ld) in %ld tapebuffers\n", reent_total_samples, total_samples, tapebuffer.size());
        }
    }
    // end take care of any samples left over in the guppi buffer

    // read blocks until we run out of data or the numblocks request is met
    int round = 0;
    do{    
        //if(control_block.vflag>=1) fprintf(stderr, "ROUND %d ========================================================\n", round);
        log_messages.printf(SCHED_MSG_LOG::MSG_DEBUG,"ROUND %d ========================================================\n", round); 
        // retval = 1 : we have a new block, which may be empty if we are fast forwarding to the point of resumption
        // retval = 0 : we do not have a new block. 
        // rawinput.invalid = 1 : we are out of raw data or have detected corruption
        // rawinput.invalid = 0 : good to continue with current file or open the next file in the series.
        retval = read_block(
	                control_block,
                    rawinput,
                    header_buf,
                    filname,
                    start_data_time,
                    end_data_time,
                    filecnt,
                    this_tapebuffer
                 );
            if(retval && this_tapebuffer.data.size() > 0) {         // if retval (OK) and no data, assume fast forward resumption
                this_tapebuffer.header.data_size = last_num_sample_bytes_read;
                tapebuffer.push_back(this_tapebuffer);
                if(control_block.vflag>=1) print_header(this_tapebuffer);
                round++;
                this_tapebuffer.data.clear();       
            }
        if(control_block.vflag>=2) {
            reent_total_samples = 0;
        }          
    } while((!(rawinput.invalid))                                                       && 
            total_samples < limit_samples);
    fprintf(stderr, "========================================================\n");
    // end read blocks until we run out of data or the numblocks request is met

    sample_bytes_left = control_block.chanbytes - last_num_sample_bytes_read;             // set up for next call

    // optionally dump a full tapebuffer
    if(dumpraw) dump_tapebuffer(file_prefix, channel);
         
    if(control_block.vflag>=1) {
             fprintf(stderr, "NUM SAMPLES total %ld final buffer contains %ld of a possible %d leaving %ld for the next time through\n", 
             total_samples, last_num_sample_bytes_read, control_block.chanbytes, control_block.chanbytes-last_num_sample_bytes_read);
    }

    if(control_block.vflag>=2) {
        reent_total_samples = 0;
        for(int buffer_i=0; buffer_i < tapebuffer.size(); buffer_i++) {
            reent_total_samples += tapebuffer[buffer_i].data.size();
            fprintf(stderr, "finish tapebuffer %d has %ld samples\n", buffer_i, tapebuffer[buffer_i].data.size());
        }
        fprintf(stderr, "finish - total_samples = %ld (%ld) in %ld tapebuffers\n", reent_total_samples, total_samples, tapebuffer.size());
    }

    fprintf(stderr, "getting coordinate history...\n");
    get_coord_history();
	
    fprintf(stderr, "finishing up...\n");

    int blocks_read = tapebuffer.size();
    if(control_block.vflag>=1) fprintf(stderr, "grabbed %ld samples from %Ld blocks covering %lf seconds of time.  Current block is %ld = %d\n", 
                                       total_samples, blocks_read, (end_data_time-start_data_time)*86400, 
                                       control_block.currentblock, tapebuffer[tapebuffer.size()-1].header.dataseq);
    start_data_time = end_data_time;    // set up for next call

    this_tapebuffer.data.clear();       
    if(total_samples != limit_samples) {
        return 0;                           // error or EOF TODO should differentiate
    } else {
        return total_samples;
    }

}   // end read_blocks_guppi()

//-------------------------------------------------------
unsigned long unpack_samples_2bit(unsigned char * raw, long int count, int pol, dr2_compact_block_t &this_tapebuffer) {
//-------------------------------------------------------
// unpack guppi 2 bit complex samples into a vector of complex signed chars.  Here we have 1 byte per sample.
// pol is polarization (0 for X, 1 for Y) 

	 long i;
     int pol_shift = pol * 2;    // pol_shift will be 0 or 2
	 float quantlookup[4];
     const int stride = 1;       // number of bytes in each time (real,imag,2pols)
     std::complex<signed char> sample;

     // the * 10 is to keep some of the precision intact as we
     // go to signed char  
	 quantlookup[0] = 3.3358750     * 10;
	 quantlookup[1] = 1.0           * 10;      
	 quantlookup[2] = -1.0          * 10;
	 quantlookup[3] = -3.3358750    * 10;

     // short circuit in the case of getting the total number of samples that we want.
     // The raw pointer will left where it is for the next time through.
	 for(i=0; i < count && total_samples < limit_samples; i+=stride) {
          // TODO - for some data sets we may need to flip i and q.  Whether or not to do this
          // should come from a boolean in the recorder_config table.
          // real (2 bits of raw data)
		  sample.real( (signed char)quantlookup[( raw[i] >> (pol_shift * 2)      & 1)        +   // bit 0 or 4  
                                                ((raw[i] >> (pol_shift * 2 + 1)  & 1) * 2)]      // bit 1 or 5
                     ); 
		  // imag (2 bits of raw data)
		  sample.imag( (signed char)quantlookup[( raw[i] >> ((pol_shift+1) * 2)      & 1)      +   // bit 2 or 6 
                                                ((raw[i] >> ((pol_shift+1) * 2 + 1)  & 1) * 2)]    // bit 3 or 7
                     ); 
          // add this sample to our vector
          this_tapebuffer.data.push_back(sample);
          total_samples++;
    }

//fprintf(stderr, "old count : %ld    new count : %ld\n", count,  this_tapebuffer.data.size());
//fprintf(stderr, "old data  : %d %d  new data  :  %d %d\n", samples[0], samples[1], this_tapebuffer.data[0].real(), this_tapebuffer.data[0].imag());
//fprintf(stderr, "old data  : %d %d  new data  :  %d %d\n", samples[2], samples[3], this_tapebuffer.data[1].real(), this_tapebuffer.data[1].imag());
//fprintf(stderr, "old data  : %d %d  new data  :  %d %d\n", samples[count*2-2], samples[count*2-1], this_tapebuffer.data[count-1].real(), this_tapebuffer.data[count-1].imag());

	 //return i/stride;   // == i in this case
	 return i;            // return the number of bytes (not samples) read
}


//-------------------------------------------------------
unsigned long unpack_samples_8bit(unsigned char * raw, long int count, int pol, dr2_compact_block_t &this_tapebuffer) {
//-------------------------------------------------------
// unpack 8 bit data. Here we have 4 bytes per sample
// pol is polarization (0 for X, 1 for Y) 

	 long i;
     int pol_shift = pol * 2;    // pol_shift will be 0 or 2
     const int stride = 4;       // number of bytes in each time (real,imag,2pols)
     std::complex<signed char> sample;

	 for(i=0; i < count && total_samples < limit_samples; i+=stride) {
          // TODO - for some data sets we may need to flip i and q.  Whether or not to do this
          // should come from a boolean in the recorder_config table.
          sample.real(raw[i+pol_shift]);
          sample.imag(raw[i+1+pol_shift]);
          // add this sample to our vector
          this_tapebuffer.data.push_back(sample);
          total_samples++;
	 }

	 //return i/stride;
	 return i;            // return the number of bytes (not samples) read
}


//-------------------------------------------------------
unsigned long unpack_samples(unsigned char * raw, long int count, int pol, dr2_compact_block_t &this_tapebuffer) {
//-------------------------------------------------------

    if(raw_bit_depth == 2) {
        return(unpack_samples_2bit(raw, count, pol, this_tapebuffer));   	
    } else if(raw_bit_depth == 8) {
        return(unpack_samples_8bit(raw, count, pol, this_tapebuffer));   	
    } else {
	    fprintf(stderr,"ERR: unsupported bit depth, exiting...\n");
	    exit(1);
    }
}


//-------------------------------------------------------
int exists(const char *fname) {
//-------------------------------------------------------
    FILE *file;
    if ((file = fopen(fname, "r"))) {
        fclose(file);
        return 1;
    }
    return 0;
}
