/*
	Librebench -- based on the miniLZO library, Whetstone, Linpack, Blowfish.
	Copyright (C) 2020 Robin@libresecurity.com
*/

/* 	Whetstone benchmark in C.  This program is a translation of the
 	original Algol version in "A Synthetic Benchmark" by H.J. Curnow
    and B.A. Wichman in Computer Journal, Vol  19 #1, February 1976.	
*/

/*
	This file is part of the LZO real-time data compression library.

	Copyright (C) 1996-2017 Markus Franz Xaver Johannes Oberhumer
	All Rights Reserved.

	The LZO library is free software; you can redistribute it and/or
	modify it under the terms of the GNU General Public License as
	published by the Free Software Foundation; either version 2 of
	the License, or (at your option) any later version.

	The LZO library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with the LZO library; see the file COPYING.
	If not, write to the Free Software Foundation, Inc.,
	51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

	Markus F.X.J. Oberhumer
	<markus@oberhumer.com>
	http://www.oberhumer.com/opensource/lzo/
*/

/*
	Blowfish source code Copyright (c) 2008, Tom Bonner.

	Permission is hereby granted, free of charge, to any person obtaining a
	copy of this software and associated documentation files (the "Software"),
	to deal in the Software without restriction, including without limitation
	the rights to use, copy, modify, merge, publish, distribute, sublicense,
	and/or sell copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in
	all copies or substantial portions of the Software.

	Except as contained in this notice, the name(s) of the above copyright
	holders shall not be used in advertising or otherwise to promote the sale,
	use or other dealings in this Software without prior written authorisation.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
	DEALINGS IN THE SOFTWARE.
*/

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>

#ifdef _OPENMP

#include <omp.h>

#endif


#include "minilzo.h"
#include "blowfish.h"
#include "whetstoned.h"
#include "linpack.h"

#define ITERATIONS	1000000 /* 1000000 Million Whetstone instructions */

//#define POUT 1
int single_precision(void)
{

	int n1, n2, n3, n4, n6, n7, n8, n9, n10, n11;
	double	x1, x2, x3, x4, x, y, z, t, t1, t2;
	double 	e1[4];
	int		i, j, k, l; 

    struct timeval s_tv, e_tv;

	gettimeofday(&s_tv, NULL);

	/* initialize constants */
	t   =   0.499975;
	t1  =   0.50025;
	t2  =   2.0;

	/* set values of module weights */

	n1  =   0 * ITERATIONS;
	n2  =  12 * ITERATIONS;
	n3  =  14 * ITERATIONS;
	n4  = 345 * ITERATIONS;
	n6  = 210 * ITERATIONS;
	n7  =  32 * ITERATIONS;
	n8  = 899 * ITERATIONS;
	n9  = 616 * ITERATIONS;
	n10 =   0 * ITERATIONS;
	n11 =  93 * ITERATIONS;

/* MODULE 1:  simple identifiers */

	x1 =  1.0;
	x2 = x3 = x4 = -1.0;

	for(i = 1; i <= n1; i += 1) {
		x1 = ( x1 + x2 + x3 - x4 ) * t;
		x2 = ( x1 + x2 - x3 - x4 ) * t;
		x3 = ( x1 - x2 + x3 + x4 ) * t;
		x4 = (-x1 + x2 + x3 + x4 ) * t;
	}

/* MODULE 2:  array elements */

	e1[0] =  1.0;
	e1[1] = e1[2] = e1[3] = -1.0;

	for (i = 1; i <= n2; i +=1) {
		e1[0] = ( e1[0] + e1[1] + e1[2] - e1[3] ) * t;
		e1[1] = ( e1[0] + e1[1] - e1[2] + e1[3] ) * t;
		e1[2] = ( e1[0] - e1[1] + e1[2] + e1[3] ) * t;
		e1[3] = (-e1[0] + e1[1] + e1[2] + e1[3] ) * t;
	}

/* MODULE 3:  array as parameter */

	for (i = 1; i <= n3; i += 1)
	{
		register int j;

		j = 0;
		lab:
		e1[0] = (  e1[0] + e1[1] + e1[2] - e1[3] ) * t;
		e1[1] = (  e1[0] + e1[1] - e1[2] + e1[3] ) * t;
		e1[2] = (  e1[0] - e1[1] + e1[2] + e1[3] ) * t;
		e1[3] = ( -e1[0] + e1[1] + e1[2] + e1[3] ) / t2;
		j += 1;
		if (j < 6)
			goto lab;			
	}

/* MODULE 4:  conditional jumps */

	j = 1;
	for (i = 1; i <= n4; i += 1) {
		if (j == 1)
			j = 2;
		else
			j = 3;

		if (j > 2)
			j = 0;
		else
			j = 1;

		if (j < 1 )
			j = 1;
		else
			j = 0;
	}

/* MODULE 5:  omitted */

/* MODULE 6:  integer arithmetic */

	j = 1;
	k = 2;
	l = 3;

	for (i = 1; i <= n6; i += 1) {
		j = j * (k - j) * (l -k);
		k = l * k - (l - j) * k;
		l = (l - k) * (k + j);

		e1[l - 2] = j + k + l;		/* C arrays are zero based */
		e1[k - 2] = j * k * l;
	}

/* MODULE 7:  trig. functions */

	x = y = 0.5;

	for(i = 1; i <= n7; i +=1) {
		x = t * atan(t2*sin(x)*cos(x)/(cos(x+y)+cos(x-y)-1.0));
		y = t * atan(t2*sin(y)*cos(y)/(cos(x+y)+cos(x-y)-1.0));
	}

/* MODULE 8:  procedure calls */

	x = y = z = 1.0;

	for (i = 1; i <= n8; i +=1)
	{
		x  = t * (x + y);
		y  = t * (x + y);
		z = (x + y) /t2;		
	}

/* MODULE9:  array references */

	j = 1;
	k = 2;
	l = 3;

	e1[0] = 1.0;
	e1[1] = 2.0;
	e1[2] = 3.0;

	for(i = 1; i <= n9; i += 1)
	{
		e1[j] = e1[k];
		e1[k] = e1[l];
		e1[l] = e1[j];
	}

/* MODULE10:  integer arithmetic */

	j = 2;
	k = 3;

	for(i = 1; i <= n10; i +=1) {
		j = j + k;
		k = j + k;
		j = k - j;
		k = k - j - j;
	}

/* MODULE11:  standard functions */

	x = 0.75;
	for(i = 1; i <= n11; i +=1)
		x = sqrt( exp( log(x) / t1));

    gettimeofday(&e_tv, NULL);
    printf("#1 Single Precision Whetstone in [ %4.4lf ] seconds @ Loops: %d\n",(float) ((e_tv.tv_sec - s_tv.tv_sec)*1000000L + e_tv.tv_usec - s_tv.tv_usec) / 1000000L, ITERATIONS);	


	return 1;

}




/**

	@ingroup blowfish
	@defgroup blowfish_selftest Blowfish Self-Test
	@{ 

  */ 

/** @internal Test vector. */ 

typedef struct __BLOWFISH_TEST_VECTOR
{
	BLOWFISH_UCHAR	Key [ 8 ];			/*!< 8-Byte key to use in the test. */ 
	BLOWFISH_ULONG	PlainText [ 2 ];	/*!< 8-Byte block of plaintext to encipher. */ 
	BLOWFISH_ULONG	CipherText [ 2 ];	/*!< 8-Byte block of expected ciphertext. */ 

} _BLOWFISH_TEST_VECTOR;

/** @internal Part of CBC/CFB/OFB Test vector. See #_BLOWFISH_EcbTv1. */ 

static const BLOWFISH_UCHAR _BLOWFISH_Tv3Key [ ] = { 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5, 0x96, 0x87 };

/** @internal Part of CBC/CFB/OFB Test vector. See #_BLOWFISH_EcbTv1. */ 

static const BLOWFISH_ULONG _BLOWFISH_Tv3Iv [ 2 ] = { 0xfedcba98, 0x76543210 };

/** @internal Part of CBC/CFB/OFB Test vector. See #_BLOWFISH_EcbTv1. */ 

static const BLOWFISH_ULONG _BLOWFISH_Tv3PlainText [ 6 ] = { 0x37363534, 0x33323120, 0x4e6f7720, 0x69732074, 0x68652074, 0x696d6520 };

/** @internal Part of CBC/CFB/OFB Test vector. See #_BLOWFISH_EcbTv1. */ 

static const BLOWFISH_ULONG _BLOWFISH_Tv3CipherText [ ] [ 6 ] = 
{
	{ 0x6b77b4d6, 0x3006dee6, 0x05b156e2, 0x74039793, 0x58deb9e7, 0x154616d9 },	/*!< CBC mode ciphertext. */ 
	{ 0xe73214a2, 0x822139ca, 0xf26ecf6d, 0x2eb9e76e, 0x3da3de04, 0xd1517200 },	/*!< CFB mode ciphertext. */ 
	{ 0xe73214a2, 0x822139ca, 0x62b343cc, 0x5b655873, 0x10dd908d, 0x0c241b22 }	/*!< OFB mode ciphertext. */ 
};

/** @internal CBC/CFB/OFB Test modes. */ 

static const BLOWFISH_MODE _BLOWFISH_Tv3Mode [ ] = { BLOWFISH_MODE_CBC, BLOWFISH_MODE_CFB, BLOWFISH_MODE_OFB };

/** @internal Throughput test vector */ 

typedef struct __BLOWFISH_THROUGHPUT_TEST
{
	BLOWFISH_MODE	Mode;		/*!< Mode to use in the test. */ 
	BLOWFISH_UCHAR	Parallel;	/*!< Specifies whether the mode can be parallelised. */ 

} _BLOWFISH_THROUGHPUT_TEST;

/** @internal Throughput test modes */ 

static const _BLOWFISH_THROUGHPUT_TEST _BLOWFISH_ThroughputTv [ ] = { { BLOWFISH_MODE_CBC, 0x01 }, { BLOWFISH_MODE_CFB, 0x01 }, { BLOWFISH_MODE_OFB, 0x00 }, { BLOWFISH_MODE_CTR, 0x01 } };

/** @internal Specifies the duration (in seconds) for how long the throughput tests should run (must be greater than 1, and preferably a multiple of 2) */ 

#define _THROUGHPUT_DURATION			10

#define _BLOWFISH_THROUGHPUT_STREAM_LENGTH		( 128 * 1024 )

/**

	@internal

	Display function name and readable return code to stdout.

	@param FunctionName	Name of the function called.

	@param ReturnCode	Return code from the function.

	@return Return code from printf().

  */ 

static int _BLOWFISH_PrintReturnCode ( char * FunctionName, BLOWFISH_RC ReturnCode )
{
	switch ( ReturnCode )
	{
		case BLOWFISH_RC_SUCCESS:
		{
			return 0;
		}
		case BLOWFISH_RC_INVALID_PARAMETER:
		{
			return printf ( "%s()=Invalid parameter!\n", FunctionName );
		}
		case BLOWFISH_RC_INVALID_KEY:
		{
			return printf ( "%s()=Invalid key!\n", FunctionName );
		}
		case BLOWFISH_RC_WEAK_KEY:
		{
			return printf ( "%s()=Weak key!\n", FunctionName );
		}
		case BLOWFISH_RC_BAD_BUFFER_LENGTH:
		{
			return printf ( "%s()=Invalid buffer length!\n", FunctionName );
		}
		case BLOWFISH_RC_INVALID_MODE:
		{
			return printf ( "%s()=Invalid mode!\n", FunctionName );
		}
		case BLOWFISH_RC_TEST_FAILED:
		{
			return printf ( "%s()=Self-test failed!\n", FunctionName );
		}
		default:
		{
			return printf ( "%s()=Unknown error!\n", FunctionName );
		}
	}
}

/**

	@internal

	Display the name of the specified mode to stdout.

	@param Mode	Mode to display.

	@return Return code from printf().

  */ 

static int _BLOWFISH_PrintMode ( BLOWFISH_MODE Mode )
{
	switch ( Mode )
	{
		case BLOWFISH_MODE_CBC:
		{
			return printf ( "  -- Mode=Cipher block chaining (CBC)\n" );
		}
		case BLOWFISH_MODE_CFB:
		{
			return printf ( "  -- Mode=Cipher feedback (CFB)\n" );
		}
		case BLOWFISH_MODE_OFB:
		{
			return printf ( "  -- Mode=Output feedback (OFB)\n" );
		}
		case BLOWFISH_MODE_CTR:
		{
			return printf ( "  -- Mode=Counter (CTR)\n" );
		}
		default:
		{
			return printf ( "  -- Mode=Invalid!\n" );
		}
	}
}

/**

	@internal

	Display a buffer as hex to stdout.

	@param Name		Name of the buffer.

	@param Buffer	Pointer to a buffer of data to display.

	@param Buffer	Length of the buffer.

  */ 

static void _BLOWFISH_PrintBuffer ( char * Name, BLOWFISH_PUCHAR Buffer, BLOWFISH_ULONG BufferLength )
{
	BLOWFISH_ULONG	i;

	printf ( "%s=0x", Name );

	for ( i = 0; i < BufferLength; i++ )
	{
		printf ( "%02x", Buffer [ i ] );
	}

	printf ( " (%d bytes)\n", (int)BufferLength );

	return;
}


static BLOWFISH_RC _BLOWFISH_Test_CBC_CFB_OFB ( BLOWFISH_PUCHAR Key, BLOWFISH_ULONG KeyLength, BLOWFISH_MODE Mode, BLOWFISH_PUCHAR PlainTextBuffer, BLOWFISH_PUCHAR CipherTextBuffer, BLOWFISH_ULONG BufferLength )
{
	BLOWFISH_RC			ReturnCode = BLOWFISH_RC_SUCCESS;
	BLOWFISH_CONTEXT	Context;
	BLOWFISH_PUCHAR		Buffer = 0;

	/* Initialise blowfish */ 

	ReturnCode = BLOWFISH_Init ( &Context, Key, KeyLength, Mode, _BLOWFISH_Tv3Iv [ 0 ], _BLOWFISH_Tv3Iv [ 1 ]  );

	_BLOWFISH_PrintReturnCode ( "BLOWFISH_Init", ReturnCode );

	/* Print key information */ 

	//_BLOWFISH_PrintMode ( Mode );

	if ( ReturnCode == BLOWFISH_RC_SUCCESS )
	{
		Buffer = (BLOWFISH_PUCHAR)malloc ( BufferLength );

		if ( Buffer != 0 )
		{
			/* Encipher the plaintext buffer */ 

			ReturnCode = BLOWFISH_EncipherBuffer ( &Context, PlainTextBuffer, Buffer, BufferLength );

			_BLOWFISH_PrintReturnCode ( "BLOWFISH_EncipherBuffer", ReturnCode );


			if ( ReturnCode == BLOWFISH_RC_SUCCESS )
			{
				/* Is the ciphertext as expected? */ 

				if ( memcmp ( CipherTextBuffer, Buffer, BufferLength ) == 0 )
				{
					/* Decipher the ciphertext buffer */ 

					ReturnCode = BLOWFISH_DecipherBuffer ( &Context, CipherTextBuffer, Buffer, BufferLength );

					_BLOWFISH_PrintReturnCode ( "BLOWFISH_DecipherBuffer", ReturnCode );

					if ( ReturnCode == BLOWFISH_RC_SUCCESS )
					{
						/* Is the plaintext as expected? */ 

						if ( memcmp ( PlainTextBuffer, Buffer, BufferLength ) != 0 )
						{
							/* Failed to decipher properly */ 

							_BLOWFISH_PrintBuffer ( "Invalid plaintext", Buffer, BufferLength );

							ReturnCode = BLOWFISH_RC_TEST_FAILED;
						}
					}
				}
				else
				{
					/* Failed to encipher properly */ 

					_BLOWFISH_PrintBuffer ( "Ciphertext", CipherTextBuffer, BufferLength );

					_BLOWFISH_PrintBuffer ( "Invalid ciphertext", Buffer, BufferLength );

					ReturnCode = BLOWFISH_RC_TEST_FAILED;
				}
			}

			free ( Buffer );
		}
	}

	/* Overwrite the blowfish context record */ 

	BLOWFISH_Exit ( &Context );

	return ReturnCode;
}

/**

	@internal

	Perform a stream based throughput test for the selected mode.

	@param Mode				Mode to use for the test.

	@param StreamBlockSize	Size of each chunk of the stream.

	@remarks Enciphers/deciphers data as a stream for #_THROUGHPUT_DURATION seconds, and then calculates the throughput.

	@return #BLOWFISH_RC_SUCCESS	Test passed successfully.

	@return Specific return code, see #BLOWFISH_RC.

  */ 

static BLOWFISH_RC _BLOWFISH_Test_Throughput ( BLOWFISH_MODE Mode, BLOWFISH_ULONG StreamBlockSize )
{
	BLOWFISH_RC			ReturnCode = BLOWFISH_RC_SUCCESS;
	BLOWFISH_CONTEXT	EncipherContext;
	BLOWFISH_CONTEXT	DecipherContext;
	BLOWFISH_PUCHAR		PlainTextBuffer = 0;
	BLOWFISH_PUCHAR		CipherTextBuffer = 0;
	BLOWFISH_SIZE_T		i = 0;
	BLOWFISH_SIZE_T		j = 0;
	BLOWFISH_ULONG		Sum = 0;
	clock_t				StartTime = 0;
	clock_t				EndTime = 0;
	clock_t				ElapsedEncipherTime = 0;
	clock_t				ElapsedDecipherTime = 0;
	float				BlocksProcessed = 0;
#ifdef _OPENMP
	BLOWFISH_SIZE_T		Threads = omp_get_max_threads ( );
#else
	BLOWFISH_SIZE_T		Threads = 1;
#endif

	/* Initialise blowfish for the encipher stream */ 

	ReturnCode = BLOWFISH_Init ( &EncipherContext, (BLOWFISH_PUCHAR)"0123456789abcdef", 16, Mode, _BLOWFISH_Tv3Iv [ 0 ], _BLOWFISH_Tv3Iv [ 1 ] );

	_BLOWFISH_PrintReturnCode ( "BLOWFISH_Init", ReturnCode );

	_BLOWFISH_PrintMode ( Mode );

	if ( ReturnCode == BLOWFISH_RC_SUCCESS )
	{
		/* Clone the context for the decipher stream */ 

		ReturnCode = BLOWFISH_CloneContext ( &EncipherContext, &DecipherContext );

		_BLOWFISH_PrintReturnCode ( "BLOWFISH_CloneContext", ReturnCode );

		if ( ReturnCode == BLOWFISH_RC_SUCCESS )
		{
			/* Allocate the plaintext buffer */ 

			PlainTextBuffer = (BLOWFISH_PUCHAR)malloc ( StreamBlockSize );

			if ( PlainTextBuffer != 0 )
			{
				/* Allocate the ciphertext buffer */ 

				CipherTextBuffer = (BLOWFISH_PUCHAR)malloc ( StreamBlockSize );

				if ( CipherTextBuffer != 0 )
				{
					/* Create original plaintext buffer (use 0x01 for ease of vectorised verification) */ 

					memset ( PlainTextBuffer, 0x01, StreamBlockSize );

					/* Begin the stream for enciphering */ 

					ReturnCode = BLOWFISH_BeginStream ( &EncipherContext );

					_BLOWFISH_PrintReturnCode ( "BLOWFISH_BeginStream", ReturnCode );

					if ( ReturnCode == BLOWFISH_RC_SUCCESS )
					{
						/* Begin the stream for deciphering */ 

						ReturnCode = BLOWFISH_BeginStream ( &DecipherContext );

						_BLOWFISH_PrintReturnCode ( "BLOWFISH_BeginStream", ReturnCode );

						if ( ReturnCode == BLOWFISH_RC_SUCCESS )
						{
							/* While the total elapsed time (in seconds) has not passed the duration threshold, keep enciphering/deciphering the stream */ 

							for ( i = 0; ( ( ElapsedEncipherTime + ElapsedDecipherTime ) / Threads ) / CLOCKS_PER_SEC <= _THROUGHPUT_DURATION; i++ )
							{
								/* Clear the ciphertext buffer */ 

								memset ( CipherTextBuffer, 0x00, StreamBlockSize );

								/* Encipher the plaintext */ 

								StartTime = clock ( );

								ReturnCode = BLOWFISH_EncipherStream ( &EncipherContext, PlainTextBuffer, CipherTextBuffer, StreamBlockSize );

								EndTime = clock ( );

								/* Compute the time elapsed enciphering (in milliseconds) */ 

								ElapsedEncipherTime += ( EndTime - StartTime );

								if ( ReturnCode != BLOWFISH_RC_SUCCESS )
								{
									_BLOWFISH_PrintReturnCode ( "BLOWFISH_EncipherStream", ReturnCode );

									break;
								}

								/* Clear the plaintext buffer */ 

								memset ( PlainTextBuffer, 0x00, StreamBlockSize );

								/* Decipher the ciphertext */ 

								StartTime = clock ( );

								ReturnCode = BLOWFISH_DecipherStream ( &DecipherContext, CipherTextBuffer, PlainTextBuffer, StreamBlockSize );

								EndTime = clock ( );

								/* Compute the time elapsed deciphering (in milliseconds) */ 

								ElapsedDecipherTime += ( EndTime - StartTime );

								if ( ReturnCode != BLOWFISH_RC_SUCCESS )
								{
									_BLOWFISH_PrintReturnCode ( "BLOWFISH_DecipherStream", ReturnCode );

									break;
								}

								/* Verify the integrity of the deciphered plaintext (sum of all bytes should equal the buffer length) */ 

								for ( j = 0, Sum = 0; j < (BLOWFISH_SIZE_T)StreamBlockSize; j++ )
								{
									Sum += PlainTextBuffer [ j ];
								}

								if ( Sum != StreamBlockSize )
								{
									printf ( "BLOWFISH_EncipherStream()/BLOWFISH_DecipherStream()=Integrity check failed!\n" );

									ReturnCode = BLOWFISH_RC_TEST_FAILED;

									break;
								}
							}

							if ( ReturnCode == BLOWFISH_RC_SUCCESS )
							{
								/* Ensure we managed to encipher/decipher enough data to determine the throughput in MB/s */ 

								if ( ( i * StreamBlockSize ) > ( 1024 * 1024 ) )
								{
									/* Adjust elapsed time based on thread count */

									ElapsedEncipherTime = ElapsedEncipherTime / Threads;
									ElapsedDecipherTime = ElapsedDecipherTime / Threads;

									/* Convert elapsed time from milliseconds to seconds (minimum 1 second) */ 

									ElapsedEncipherTime = ElapsedEncipherTime / CLOCKS_PER_SEC > 1 ? ElapsedEncipherTime / CLOCKS_PER_SEC : 1;
									ElapsedDecipherTime = ElapsedDecipherTime / CLOCKS_PER_SEC > 1 ? ElapsedDecipherTime / CLOCKS_PER_SEC : 1;

									/* Calculate and display the stream length */ 

									BlocksProcessed = (float)( (float)( i * StreamBlockSize ) / ( 1024 * 1024 ) );

									printf ( "  -- Stream length=%0.2f MB (%d*%d byte blocks)\n", BlocksProcessed, (int)i, (int)StreamBlockSize );

									/* Calculate and display throughputs */ 

									printf ( "  -- Encipher throughput=%0.2f MB/s\n", BlocksProcessed / ElapsedEncipherTime );
									printf ( "  -- Decipher throughput=%0.2f MB/s\n", BlocksProcessed / ElapsedDecipherTime );
								}
								else
								{
									printf ( "Failed to process enough data to determine throughput in MB/s!\n" );

									ReturnCode = BLOWFISH_RC_TEST_FAILED;
								}
							}

							/* Finished decipering */ 

							BLOWFISH_EndStream ( &DecipherContext );
						}

						/* Finished encipering */ 

						BLOWFISH_EndStream ( &EncipherContext );
					}

					/* Free the ciphertext */ 

					free ( CipherTextBuffer );
				}

				/* Free the plaintext */ 

				free ( PlainTextBuffer );
			}
		}

		/* Overwrite the deciper stream context record */ 

		BLOWFISH_Exit ( &DecipherContext );
	}

	/* Overwrite the encipher stream context record */ 

	BLOWFISH_Exit ( &EncipherContext );

	return ReturnCode;
}

/**

	@internal

	Perform all self-tests for all supported modes.

	@return #BLOWFISH_RC_SUCCESS	All tests passed successfully.

	@return Specific return code, see #BLOWFISH_RC.

  */ 

static BLOWFISH_RC _BLOWFISH_SelfTest ( )
{
	BLOWFISH_RC		ReturnCode;
	BLOWFISH_ULONG	i = 0;
	struct timeval s_tv, e_tv;

	/* Perform CBC, CFB and OFB tests on test vector 3 */ 

	for ( i = 0; i < sizeof ( _BLOWFISH_Tv3Mode ) / sizeof ( _BLOWFISH_Tv3Mode [ 0 ] ); i++ )
	{
		ReturnCode = _BLOWFISH_Test_CBC_CFB_OFB ( (BLOWFISH_PUCHAR)&_BLOWFISH_Tv3Key, sizeof ( _BLOWFISH_Tv3Key ), _BLOWFISH_Tv3Mode [ i ], (BLOWFISH_PUCHAR)&_BLOWFISH_Tv3PlainText, (BLOWFISH_PUCHAR)&_BLOWFISH_Tv3CipherText [ i ], sizeof ( _BLOWFISH_Tv3PlainText ) );

		if ( ReturnCode != BLOWFISH_RC_SUCCESS )
		{
			return ReturnCode;
		}
	}

#ifdef _OPENMP

	gettimeofday(&s_tv, NULL);
	/* Perform parallelised throughput tests if there is more than 1 available thread */ 

	if ( omp_get_max_threads ( ) > 1 )
	{
		printf ( "#6 Parallelised blowfish bench (using %d threads for ~%d seconds per mode)\n", omp_get_max_threads ( ), _THROUGHPUT_DURATION );

		for ( i = 0; i < sizeof ( _BLOWFISH_ThroughputTv ) / sizeof ( _BLOWFISH_ThroughputTv [ 0 ] ); i++ )
		{
			if ( _BLOWFISH_ThroughputTv [ i ].Parallel != 0x00 )
			{
				ReturnCode = _BLOWFISH_Test_Throughput ( _BLOWFISH_ThroughputTv [ i ].Mode, _BLOWFISH_THROUGHPUT_STREAM_LENGTH );

				if ( ReturnCode != BLOWFISH_RC_SUCCESS )
				{
					return ReturnCode;
				}
			}
		}
	}
	/* Use a single thread for the serialised tests */ 
	gettimeofday(&e_tv, NULL);
	printf("  -- Parallelised blowfish bench in [ %4.4lf ] seconds\n", (float) ((e_tv.tv_sec - s_tv.tv_sec)*1000000L + e_tv.tv_usec - s_tv.tv_usec) / 1000000L);

	omp_set_num_threads ( 1 );

#endif

	/* Perform serialised throughput tests */ 

	printf ( "#7 Serialised blowfish bench (using 1 thread for ~%d seconds per mode)\n", _THROUGHPUT_DURATION );
	gettimeofday(&s_tv, NULL);

	for ( i = 0; i < sizeof ( _BLOWFISH_ThroughputTv ) / sizeof ( _BLOWFISH_ThroughputTv [ 0 ] ); i++ )
	{
		ReturnCode = _BLOWFISH_Test_Throughput ( _BLOWFISH_ThroughputTv [ i ].Mode, _BLOWFISH_THROUGHPUT_STREAM_LENGTH );

		if ( ReturnCode != BLOWFISH_RC_SUCCESS )
		{
			return ReturnCode;
		}
	}
	gettimeofday(&e_tv, NULL);
	printf("  -- Serialised blowfish bench in [ %4.4lf ] seconds\n", (float) ((e_tv.tv_sec - s_tv.tv_sec)*1000000L + e_tv.tv_usec - s_tv.tv_usec) / 1000000L);
	return BLOWFISH_RC_SUCCESS;
}

/**

	@internal

	Main entry point for the blowfish self-test application.

	@param ArgumentCount	Number of command line arguments passed to the application.

	@param ArgumentVector	Array of command line arguments passed to the application.

	@return #BLOWFISH_RC_SUCCESS	All tests passed successfully.

	@return Specific return code, see #BLOWFISH_RC.

  */ 

/* We want to compress the data block at 'in' with length 'IN_LEN' to
 * the block at 'out'. Because the input block may be incompressible,
 * we must provide a little more output space in case that compression
 * is not possible.
 */

#define IN_LEN      (1280*1024000ul)
#define OUT_LEN     (IN_LEN + IN_LEN / 16 + 64 + 3)

static unsigned char __LZO_MMODEL in  [ IN_LEN ];
static unsigned char __LZO_MMODEL out [ OUT_LEN ];


/* Work-memory needed for compression. Allocate memory in units
 * of 'lzo_align_t' (instead of 'char') to make sure it is properly aligned.
 */

#define HEAP_ALLOC(var,size) \
    lzo_align_t __LZO_MMODEL var [ ((size) + (sizeof(lzo_align_t) - 1)) / sizeof(lzo_align_t) ]

static HEAP_ALLOC(wrkmem, LZO1X_1_MEM_COMPRESS);

void lzo_bench(void)
{
	int r;
    lzo_uint in_len;
    lzo_uint out_len;
    lzo_uint new_len;
	struct timeval s_tv, e_tv;
/*
 * Step 1: initialize the LZO library
 */
    if (lzo_init() != LZO_E_OK)
    {
        printf("internal error - lzo_init() failed !!!\n");
    }

/*
 * Step 2: prepare the input block that will get compressed.
 *         We just fill it with zeros in this example program,
 *         but you would use your real-world data here.
 */
    in_len = IN_LEN;
    lzo_memset(in,rand() % 256,in_len);
    
    gettimeofday(&s_tv, NULL);
/*
 * Step 3: compress from 'in' to 'out' with LZO1X-1
 */
    r = lzo1x_1_compress(in,in_len,out,&out_len,wrkmem);
	r = lzo1x_1_compress(in,in_len,out,&out_len,wrkmem);
    r = lzo1x_1_compress(in,in_len,out,&out_len,wrkmem);
    r = lzo1x_1_compress(in,in_len,out,&out_len,wrkmem);
    r = lzo1x_1_compress(in,in_len,out,&out_len,wrkmem);			
    
    if(r != LZO_E_OK)
    {
        /* this should NEVER happen */
        printf("internal error - compression failed: %d\n", r);
    }
    /* check for an incompressible block */
    if (out_len >= in_len)
    {
        printf("This block contains incompressible data.\n");
    }

/*
 * Step 4: decompress again, now going from 'out' to 'in'
 */
    new_len = in_len;
    r = lzo1x_decompress(out,out_len,in,&new_len,NULL);
    r = lzo1x_decompress(out,out_len,in,&new_len,NULL);
    r = lzo1x_decompress(out,out_len,in,&new_len,NULL);
	r = lzo1x_decompress(out,out_len,in,&new_len,NULL);		
    r = lzo1x_decompress(out,out_len,in,&new_len,NULL);	
    
    if (r == LZO_E_OK && new_len == in_len)
    {
        gettimeofday(&e_tv, NULL);
        printf("#4 Decompressed %lu MB to %lu MB, Lzo bench in [ %4.4lf ] seconds\n",(unsigned long) out_len/(1024*1024), 
            (unsigned long) in_len/(1024*1024), (float) ((e_tv.tv_sec - s_tv.tv_sec)*1000000L + e_tv.tv_usec - s_tv.tv_usec) / 1000000L);
    }       
    else
    {
        printf("internal error - decompression failed: %d\n", r);
    }
}

void blowfish_bench(void)
{
	_BLOWFISH_SelfTest ( );

}

void superpi_bench(void)
{
    long loops = 3000000000L;
	double a = 1.0f;
	double b = 1.0f / sqrt(2);
	double t = 1.0f / 4.0f;
	double p = 1.0f;
	struct timeval s_tv, e_tv;

	gettimeofday(&s_tv, NULL);
	for (long i = 0; i < loops; i++)
	{
		double x = (a + b) / 2;
		double y = sqrt(a * b);
		double w = t - p * pow(a - x, 2);
		double z = 2 * p;
		a = x;
		b = y;
		t = w;
		p = z;
	}	
	gettimeofday(&e_tv, NULL);
	printf("#5 Superpi bench in [ %4.4lf ] seconds\n", (float) ((e_tv.tv_sec - s_tv.tv_sec)*1000000L + e_tv.tv_usec - s_tv.tv_usec) / 1000000L);

}

int main(void)
{

	printf("-- Libre Bench --\n");

	single_precision();
    
	double_precision();
    
	linpack_bench();

	lzo_bench();

	superpi_bench();
	
	blowfish_bench();

   
    return 0;
}
