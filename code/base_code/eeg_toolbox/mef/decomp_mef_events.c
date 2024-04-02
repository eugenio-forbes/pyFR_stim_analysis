//
//mex decomp_mef_events.c RED_decode.c RED_encode.c endian_functions.c mef_lib.c AES_encryption.c crc_32.c

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "mef_header.h"
#include "mef_lib.h"
#include "RED_codec.h"
#include "size_types.h"

#define BIG_ENDIAN_CODE		0
#define LITTLE_ENDIAN_CODE	1

void decomp_mef_events(char *f_name, ui4 n_events, ui8 *start_idx_array, 
                       ui4 duration, char *password, si4 *decomp_data_all, ui4 *decomp_data_lengths )
{
	char    *c, *comp_data, *cdp, *last_block_p, encryptionKey[240];
	ui1		*header, *diff_buffer;
    ui1     cpu_endianness();
	si4     *decomp_data, *dcdp, *temp_data_buf;
	ui4		n_read, comp_data_len, bytes_decoded, tot_samples;
	ui4		i, iE, n_index_entries, tot_index_fields, kept_samples, skipped_samples;
	ui8     start_block_file_offset, end_block_file_offset, start_block_idx, end_block_idx;
	ui8     *index_data, last_block_len, RED_decompress_block();
    ui8     start_idx, end_idx;
    
	FILE                *fp;
	MEF_HEADER_INFO		hdr_info;
	RED_BLOCK_HDR_INFO	block_hdr;
    
	void AES_KeyExpansion();
		
	// get cpu endianness
	if (cpu_endianness() != LITTLE_ENDIAN_CODE) {
		mexErrMsgTxt("[decomp_mef_events] is currently only compatible with little-endian machines => exiting");
		return;
	}
	
	// read header
	fp = fopen(f_name, "r");
	if (fp == NULL) { 
		printf("[decomp_mef_events] could not open the file \"%s\" => exiting\n",  f_name);
		return;
	}
	header = (ui1 *) malloc(MEF_HEADER_LENGTH);  // malloc to ensure boundary alignment
	n_read = fread((void *) header, sizeof(ui1), (size_t) MEF_HEADER_LENGTH, fp);
	if (n_read != MEF_HEADER_LENGTH) {
		printf("[decomp_mef_events] error reading the file \"%s\" => exiting\n",  f_name);
		return;
	}	
	if ((read_mef_header_block(header, &hdr_info, password))) {
		printf("[decomp_mef_events] header read error for file \"%s\" => exiting\n", f_name);
		return;		
	}
	free(header);
	
	// showHeader(&hdr_info);
	
	// get file endianness
	if (hdr_info.byte_order_code != LITTLE_ENDIAN_CODE) {
		mexErrMsgTxt("[decomp_mef_events] is currently only compatible with little-endian files (file \"%s\") => exiting");
		return;
	}

	if (hdr_info.data_encryption_used) {
		AES_KeyExpansion(4, 10, encryptionKey, hdr_info.session_password); 
	}
	else
		*encryptionKey = 0;
	
	// read in index data
	n_index_entries = (ui4) hdr_info.number_of_index_entries;
	fseeko(fp, (off_t) hdr_info.index_data_offset, SEEK_SET);
	tot_index_fields = n_index_entries * 3;	// 3 fields per entry
	index_data = (ui8 *) malloc(tot_index_fields * sizeof(ui8));
	if (index_data == NULL) {
		printf("[decomp_mef_events] could not allocate enough memory for file \"%s\" => exiting\n", f_name);
		return;
	}
	
	n_read = fread(index_data, sizeof(ui8), (size_t) tot_index_fields, fp);
	if (n_read != tot_index_fields) {
		printf("[decomp_mef_events] error reading index data for file \"%s\" => exiting\n", f_name);
		return;
	}
    
    // allocate input buffer
    comp_data_len = (ui4) ((duration + 2*hdr_info.maximum_block_length) * sizeof(si4));
    comp_data = (char *) malloc(comp_data_len);
    if (comp_data == NULL) {
        printf("[decomp_mef_events] could not allocate enough memory for file \"%s\" => exiting\n", f_name);
        return;
    }
    
    diff_buffer = (ui1 *) malloc(hdr_info.maximum_block_length * sizeof(si4));
    temp_data_buf = (si4 *) malloc(hdr_info.maximum_block_length * sizeof(si4));
    
    // loop over all events
    for (iE = 0; iE < n_events; iE++) {
        start_idx = start_idx_array[iE];
        end_idx = start_idx + duration - 1;
        decomp_data = decomp_data_all + iE*duration;

        // find block containing start of requested range
        if (start_idx >= hdr_info.number_of_samples) {
//             printf("[decomp_mef_events] start index for file \"%s\" exceeds the number of samples in the file => exiting\n", f_name);
            continue;
        }
        for (i = 2; i < tot_index_fields; i += 3) {
            //printf("i %d\tidx data %lu\n", i, index_data[i-1]);
            if (index_data[i] > start_idx)
                break;
        }
        i -= 3; // rewind one triplet
        start_block_idx = index_data[i]; // sample index of start of block containing start index
        start_block_file_offset = index_data[i - 1];  // file offset of block containing start index

        // find block containing end of requested range
        if (end_idx >= hdr_info.number_of_samples) {
//             printf("[decomp_mef_events] end index for file \"%s\" exceeds the number of samples in the file => tail values will be zeros\n", f_name);
            end_idx = hdr_info.number_of_samples - 1;
        }
        for (; i < tot_index_fields; i += 3)
            if (index_data[i] > end_idx)
                break;
        i -= 3; // rewind one triplet
        end_block_idx = index_data[i]; // sample index of start of block containing end index
        end_block_file_offset = index_data[i - 1];  // file offset of block containing end index
        if (i == (tot_index_fields - 1))
            last_block_len = hdr_info.index_data_offset - end_block_file_offset;  // file offset of index data
        else
            last_block_len = index_data[i + 2] - end_block_file_offset;  // file offset of next block

        // read in compressed data
        comp_data_len = (ui4) (end_block_file_offset - start_block_file_offset + last_block_len);
        fseeko(fp, (off_t) start_block_file_offset, SEEK_SET);
        n_read = fread(comp_data, sizeof(char), (size_t) comp_data_len, fp);
        if (n_read != comp_data_len) {
            printf("[decomp_mef_events] error reading data for file \"%s\" => exiting\n", f_name);
            break;
        }

        // decompress data

        // decode first block to temp array
        cdp = comp_data;  
        bytes_decoded = (ui4) RED_decompress_block(cdp, temp_data_buf, diff_buffer, encryptionKey, 0, &block_hdr);
        cdp += bytes_decoded;

        // copy requested samples from first block to output buffer
        skipped_samples = (ui4) (start_idx - start_block_idx);
        if (skipped_samples > block_hdr.sample_count) {
            //this is bad- likely means idx data is corrupt 
            printf("[decomp_mef_events] block indexing error: decoded %d samples, attempting to skip %lu samples\n", 
                    block_hdr.sample_count, skipped_samples);
            continue;
        }
        
        kept_samples = block_hdr.sample_count - skipped_samples;
        tot_samples = (ui4) (end_idx - start_idx + 1);
        if (kept_samples >= tot_samples)
        {
            // start and end indices in same block => already done
            memcpy((void *) decomp_data, (void *) (temp_data_buf + skipped_samples), tot_samples * sizeof(si4));
        }
        else 
        {
            memcpy((void *) decomp_data, (void *) (temp_data_buf + skipped_samples), kept_samples * sizeof(si4));
            tot_samples = kept_samples;

            dcdp = decomp_data + kept_samples;
            last_block_p = comp_data + (ui4) (end_block_file_offset - start_block_file_offset);
            while (cdp < last_block_p) {
                bytes_decoded = (ui4) RED_decompress_block(cdp, dcdp, diff_buffer, encryptionKey, 0, &block_hdr);
                cdp += bytes_decoded;
                dcdp += block_hdr.sample_count; 
            }

            // decode last block to temp array
            (void) RED_decompress_block(cdp, temp_data_buf, diff_buffer, encryptionKey, 0, &block_hdr);

            // copy requested samples from last block to output buffer
            kept_samples = (ui4) (end_idx - end_block_idx + 1);
            memcpy((void *) dcdp, (void *) temp_data_buf, kept_samples * sizeof(si4));
            tot_samples += kept_samples;
        }
        
        decomp_data_lengths[iE] = tot_samples;
    }
    
    fclose(fp);
    
    free(index_data);
	free(comp_data);
	free(diff_buffer);
	free(temp_data_buf);

	return;
}


// The mex gateway routine 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char    *f_name, *password;
	si4     buf_len, status, dims[2], *eeg_data;
    ui4     n_events, duration, *eeg_data_lengths;
    sf8     *start_idx_arrayF;
	ui8     *start_idx_array;
    int     i;
	void	decomp_mef_events();
	
	//  Check for proper number of arguments 
	if ( nrhs < 4 || nrhs > 4 || nlhs < 1 || nlhs > 2 ) 
        mexErrMsgTxt("[decomp_mef_events] invalid number of input or output arguments. Type 'help decomp_mef_events' for help.");
	
	// get the input file name (argument 1)
	if (mxIsChar(prhs[0]) != 1) { // Check to make sure the first input argument is a string 
		mexErrMsgTxt("[decomp_mef_events] file name must be a string => exiting");
		return;
	}	
	buf_len = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 2; // Get the length of the input string 
	f_name = malloc(buf_len); // Allocate memory for file_name string
	status = mxGetString(prhs[0], f_name, buf_len);
	if (status != 0) {
		mexWarnMsgTxt("[decomp_mef_events] not enough space for input file name string => exiting");
		return;
	}
    
    //  get the start index array (argument 2)
    n_events = (ui4) (mxGetN(prhs[1]) * mxGetM(prhs[1]));
    if( mxIsUint64(prhs[1]) ) {
        start_idx_array = (ui8*) mxGetData(prhs[1]);
    } else if( mxIsDouble(prhs[1]) ) {
        start_idx_arrayF = (sf8*) mxGetPr(prhs[1]);
        start_idx_array = (ui8*) malloc(n_events * sizeof(ui8));
        for( i = 0; i < n_events; i++ )
            start_idx_array[i] = (ui8) start_idx_arrayF[i];
    } else {
        mexErrMsgTxt("[decomp_mef_events] start index array must be in double or uint64 format => exiting");
    }
	
	//  get the duration (argument 3)
	if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetN(prhs[2]) * mxGetM(prhs[2]) != 1) ) {
		mexErrMsgTxt("[decomp_mef_events] duration must be a scalar => exiting");
		return;
	}	
	duration = (ui4) mxGetScalar(prhs[2]);

	// get the password (argument 4)
	if (mxIsChar(prhs[3]) != 1) { // Check to make sure the fourth input argument is a string 
		mexErrMsgTxt("[decomp_mef_events] Password must be a stringx => exiting");
		return;
	}	
	buf_len = (mxGetM(prhs[3]) * mxGetN(prhs[3])) + 2; // Get the length of the input string 
	password = malloc(buf_len); // Allocate memory for file_name string
	status = mxGetString(prhs[3], password, buf_len);
	if (status != 0) {
		mexWarnMsgTxt("[decomp_mef_events] not enough space for password string => exiting");
		return;
	}

	// Set the output pointer to the eeg_data output matrix. 
	dims[0] = duration; dims[1] = n_events;
	plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
	
	// Create a C pointer to a copy of the output matrix. 
	eeg_data = (si4 *) mxGetData(plhs[0]);
	if (eeg_data == NULL) {
		mexErrMsgTxt("[decomp_mef_events] could not allocate enough memory => exiting");
		return;
	}
    
    // Set the output pointer to the output matrix. 
    dims[0] = n_events; dims[1] = 1;
	plhs[1] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
	
	// Create a C pointer to a copy of the output matrix. 
	eeg_data_lengths = (ui4 *) mxGetData(plhs[1]);
	if (eeg_data_lengths == NULL) {
		mexErrMsgTxt("[decomp_mef_events] could not allocate enough memory => exiting");
		return;
	}
	
	// Call the C subroutine. 
	decomp_mef_events(f_name, n_events, start_idx_array, duration, password, eeg_data, eeg_data_lengths);
	
	return;
} 
