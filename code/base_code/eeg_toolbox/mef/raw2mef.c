/*

 *  raw2mef.c
 *  
 
 Multiscale electrophysiology format example program
 Program to read raw 32bit integers and save data in mef format (v.2) 
 
 To compile for a 64-bit intel system: (options will vary depending on your particular compiler and platform)
 Intel Compiler: icc raw2mef.c RED_encode.c mef_lib.c endian_functions.c AES_encryption.c crc_32.c -o raw2mef -fast -m64
 GCC: gcc raw2mef.c RED_encode.c mef_lib.c endian_functions.c AES_encryption.c crc_32.c -o raw2mef -O3 -arch x86_64
 
 
 This software is made freely available under the GNU public license: http://www.gnu.org/licenses/gpl-3.0.txt
 
 Thanks to all who acknowledge the Mayo Systems Electrophysiology Laboratory, Rochester, MN
 in academic publications of their work facilitated by this software.
 
 mex raw2mef.c RED_encode.c mef_lib.c endian_functions.c AES_encryption.c crc_32.c
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#include "size_types.h"
#include "mef_header_2_0.h"
#include "RED_codec.h"

#include "mex.h"

int raw2mef( si4 *in_data, const char *outFileName, const char *subject_password, const char *session_password, MEF_HEADER_INFO *header ) {
	/*si4  *in_data, *in_data4, */
	si4 *idp, i, j, num, numBlocks, l, value, RED_block_size, num_entries, max_block_size; 
	si4 read_size, discontinuity_flag, max_value, min_value, byte_offset;
	ui8 numEntries, inDataLength, bytesDecoded, entryCounter, dataCounter, block_hdr_time, index_data_offset;
	si1 encryptionKey[240];
	ui1 *hdr_bk, *data, *dp, *diff_buffer, *dbp, byte_padding[8];
	FILE *fp;
	RED_BLOCK_HDR_INFO RED_bk_hdr;
	INDEX_DATA *index_block, *ip;
	void AES_KeyExpansion();
	unsigned long RED_compress_block();
    int status;
	
	memset(encryptionKey, 0, 240);
	/*if (*session_password) */
    if( strlen(session_password) > 0 )
	{
		/*check password length*/
		if (strlen(session_password) > 16) {
			fprintf(stderr, "Error: Password cannot exceed 16 characters\n");
			return(1);
		}
		header->session_encryption_used = 1;
		strcpy(header->session_password, session_password);
        
		/*RED block header encryption is used here: comment next two lines to disable*/
		AES_KeyExpansion(4, 10, encryptionKey, header->session_password); 
		header->data_encryption_used = 1;
        /*header->session_encryption_used = 1;*/
	}

	/*if (*subject_password)*/
    if( strlen(subject_password) > 0 )
	{
		/*check password length*/
		if (strlen(subject_password) > 16) {
			fprintf(stderr, "Error: Password cannot exceed 16 characters\n");
			return(1);
		}
		header->subject_encryption_used = 1;
	}
	
	header->maximum_block_length = header->block_interval * header->sampling_frequency;
	numBlocks = ceil((double)header->number_of_samples/header->maximum_block_length);
    
    data = malloc(4*header->number_of_samples);
    index_block = (INDEX_DATA *)calloc(3*numBlocks, sizeof(ui8));
    
    if (data == NULL || index_block == NULL) {
		fprintf(stderr, "malloc error\n");
		return 1;
	}

	dp = data;	idp = in_data; ip = index_block;
	dataCounter = 0; entryCounter = 0; discontinuity_flag = 1;
	max_value = 1<<31; min_value = max_value-1; max_block_size = 0;
	num_entries = header->maximum_block_length;
	for ( i = 0; i < numBlocks; i++ )
	{
		if (i != 0) {
			discontinuity_flag = 0; 
			ip++;
		}
		
		ip->time = header->recording_start_time + i * header->block_interval*1000000;
		ip->file_offset = dataCounter + MEF_HEADER_LENGTH; 
		ip->sample_number = i * header->maximum_block_length;
		
		if (header->number_of_samples - entryCounter < header->maximum_block_length)
			num_entries = header->number_of_samples - entryCounter;		
		
		RED_block_size = RED_compress_block(idp, dp, num_entries, ip->time, 
											(ui1)discontinuity_flag, encryptionKey, &RED_bk_hdr);
		dp += RED_block_size; 
		dataCounter += RED_block_size;
		idp += RED_bk_hdr.sample_count;
		entryCounter += RED_bk_hdr.sample_count;

		if (RED_bk_hdr.max_value > max_value) max_value = RED_bk_hdr.max_value;
		if (RED_bk_hdr.min_value < min_value) min_value = RED_bk_hdr.min_value;
		if (RED_block_size > max_block_size) max_block_size = RED_block_size;
	}
	
	/*free(in_data); in_data = NULL;*/

	header->maximum_data_value = max_value; 
	header->minimum_data_value = min_value;
	header->maximum_compressed_block_size = max_block_size;
	header->number_of_index_entries = numBlocks;
	header->index_data_offset = dataCounter + MEF_HEADER_LENGTH;
	header->recording_end_time = ip->time + (unsigned long)((double)num_entries*1000000.0/header->sampling_frequency + 0.5);
	sprintf(header->compression_algorithm, "Range Encoded Differences (RED)");
	sprintf(header->encryption_algorithm,  "AES %d-bit", ENCRYPTION_BLOCK_BITS);
	
	hdr_bk = calloc(sizeof(ui1), MEF_HEADER_LENGTH);
	memset(hdr_bk, 0, MEF_HEADER_LENGTH);
	
	fp = fopen(outFileName, "w");
	if (fp == NULL)
        mexErrMsgTxt("[raw2mef] Error opening file");

	printf("\n\nWriting file %s: %ld entries \n", outFileName, entryCounter);
	
	/*write blank header block as a place holder*/
	num = fwrite(hdr_bk, 1, MEF_HEADER_LENGTH, fp);
	
	num = fwrite(data, sizeof(ui1), dataCounter, fp);
	if (num != dataCounter) {
		fclose(fp); free(data); free(index_block);
        mexErrMsgTxt("[raw2mef] Error writing file data");
	}

	/*byte align index data if needed*/
	index_data_offset = ftell(fp);
	byte_offset = index_data_offset % 8;
	if (byte_offset) {
		memset(byte_padding, 0, 8);
		fwrite(byte_padding, sizeof(ui1), 8 - byte_offset, fp);
		index_data_offset += 8 - byte_offset;
	}
	header->index_data_offset = index_data_offset;
	
	num = fwrite(index_block, sizeof(INDEX_DATA), numBlocks, fp);
	if (num != numBlocks) {
		fclose(fp); free(data); free(index_block);
        mexErrMsgTxt("[raw2mef] Error writing file block indices");
	}
    
	fseek(fp, 0, SEEK_SET);
	status = build_mef_header_block(hdr_bk, header, subject_password);
    if (status)
        mexErrMsgTxt("[raw2mef] Error building header block");

	num = fwrite(hdr_bk, 1, MEF_HEADER_LENGTH, fp);
    if( num != MEF_HEADER_LENGTH )
        mexErrMsgTxt("[raw2mef] Error writing header block");
	
	fclose(fp);
	free(index_block);
	free(data);
	return 0;
}

/* The mex gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* inputs: inData, samplingRate, outFileName, subjectPassword, sessionPassword*/
	si4 *inData;
//     sf8 *inDataF;
	si1 *outFileName, *subjectPassword, *sessionPassword;
	MEF_HEADER_INFO header;
	
	int	buf_len, status, i;
	
    /*set defaults*/
	memset(&header, 0, sizeof(MEF_HEADER_INFO));
	header.byte_order_code = 1; /*We assume little-endian*/
	header.header_version_major = 2;
	header.header_version_minor = 0;
	header.subject_encryption_used = 0;
	header.session_encryption_used = 0;
    header.block_interval = 1; /*seconds*/
	
	/*  Check for proper number of arguments */
	if (nrhs < 5 || nrhs > 6 || nlhs != 0) 
		mexErrMsgTxt("[raw2mef] invalid number of input or output arguments. Type 'help raw2mef' for help.");
	
	/* get the input data (argument 0)*/
    header.number_of_samples = mxGetNumberOfElements( prhs[0] );
    
	if( mxIsInt32(prhs[0]) )
        inData = mxGetData( prhs[0] );
//     else if( mxIsDouble(prhs[0]) ) {
//         inDataF = (sf8*) mxGetPr(prhs[0]);
//         inData = (si4*) malloc(header.number_of_samples * sizeof(si4));
//         for( i = 0; i < header.number_of_samples; i++ )
//             inData[i] = (si4) inDataF[i];
//     } 
    else
		mexErrMsgTxt("[raw2mef] inData must be in int32 format");
	
	/* get the sampling rate (argument 1)*/
	if ( !mxIsNumeric(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1 )
		mexErrMsgTxt("[raw2mef] sampling rate must be a numeric scalar");
	
	header.sampling_frequency = (sf8) mxGetScalar( prhs[1] );

	/* get the output file name (argument 2)*/
	if ( !mxIsChar(prhs[2]) )
		mexErrMsgTxt("[raw2mef] file name must be a string");

	buf_len = mxGetNumberOfElements(prhs[2]) + 2; /* Get the length of the input string */
	outFileName = malloc(buf_len); /* Allocate memory for file_name string*/
	status = mxGetString(prhs[2], outFileName, buf_len);
	if (status)
		mexWarnMsgTxt("[raw2mef] not enough space for output file name");
		
	/* get the subject password (argument 3)*/
	if ( !mxIsChar(prhs[3]) )
		mexErrMsgTxt("[raw2mef] subject password must be a string");

	buf_len = mxGetNumberOfElements(prhs[3]) + 2; /* Get the length of the input string */
	subjectPassword = malloc(buf_len); /* Allocate memory for file_name string*/
	status = mxGetString(prhs[3], subjectPassword, buf_len);
	if (status)
		mexWarnMsgTxt("[raw2mef] not enough space for subject password");
		
	/* get the session password (argument 4)*/
	if ( !mxIsChar(prhs[4]) )
		mexErrMsgTxt("[raw2mef] session password must be a string");

	buf_len = mxGetNumberOfElements(prhs[4]) + 2; /* Get the length of the input string */
	sessionPassword = malloc(buf_len); /* Allocate memory for file_name string*/
	status = mxGetString(prhs[4], sessionPassword, buf_len);
	if (status)
		mexWarnMsgTxt("[raw2mef] not enough space for session password");
    
    /* get block interval (argument 5)*/
    if ( nrhs >= 6 ) {
        if ( !mxIsNumeric(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1 )
            mexErrMsgTxt("[raw2mef] block interval must be a numeric scalar");

        header.block_interval = (ui8) mxGetScalar( prhs[5] );
    }
	
	/* Call the C subroutine. */
	status = raw2mef( inData, outFileName, subjectPassword, sessionPassword, &header );
	
	return;
}

