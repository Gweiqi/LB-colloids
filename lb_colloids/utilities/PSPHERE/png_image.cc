#include <stdio.h>
#include <glib.h>
#include "png.h"
#include "matrix.h"

//#include "d2_quad.h"

/* This is an example of how to use libpng to read and write PNG files.
 * The file libpng.txt is much more verbose then this.  If you have not
 * read it, do so first.  This was designed to be a starting point of an
 * implementation.  This is not officially part of libpng, is hereby placed
 * in the public domain, and therefore does not require a copyright notice.
 *
 * This file does not currently compile, because it is missing certain
 * parts, like allocating memory to hold an image.  You will have to
 * supply these parts to get it to compile.  For an example of a minimal
 * working PNG reader/writer, see pngtest.c, included in this distribution;
 * see also the programs in the contrib directory.
 */

#include "png_image.h"


 /* The png_jmpbuf() macro, used in error handling, became available in
  * libpng version 1.0.6.  If you want to be able to run your code with older
  * versions of libpng, you must define the macro yourself (but only if it
  * is not already defined by libpng!).
  */
#ifndef png_jmpbuf
#  define png_jmpbuf(png_ptr) ((png_ptr)->jmpbuf)
#endif


/* Check to see if a file is a PNG file using png_sig_cmp().  png_sig_cmp()
 * returns zero if the image is a PNG and nonzero if it isn't a PNG.
 *
 * The function check_if_png() shown here, but not used, returns nonzero (true)
 * if the file can be opened and is a PNG, 0 (false) otherwise.
 *
 * If this call is successful, and you are going to keep the file open,
 * you should call png_set_sig_bytes(png_ptr, PNG_BYTES_TO_CHECK); once
 * you have created the png_ptr, so that libpng knows your application
 * has read that many bytes from the start of the file.  Make sure you
 * don't call png_set_sig_bytes() with more than 8 bytes read or give it
 * an incorrect number of bytes read, or you will either have read too
 * many bytes (your fault), or you are telling libpng to read the wrong
 * number of magic bytes (also your fault).
 *
 * Many applications already read the first 2 or 4 bytes from the start
 * of the image to determine the file type, so it would be easiest just
 * to pass the bytes to png_sig_cmp() or even skip that if you know
 * you have a PNG file, and call png_set_sig_bytes().

 Code moved to read_image()

 */
#define PNG_BYTES_TO_CHECK 4

png_image::png_image()
{
}

png_image::~png_image()
{
   /* clean up after the read, and free any memory allocated - REQUIRED */
   png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
}

guint32 png_image::read_image(char *file_name)
{
   unsigned int sig_read = PNG_BYTES_TO_CHECK;
   int png_transforms;
   int row_bytes; //,pixel_size;
   //int i;

   png_byte buf[PNG_BYTES_TO_CHECK];

   error=0;
 
   if ((fp = fopen(file_name, "rb")) == NULL)
     error=FILE_NOT_FOUND;

  /* Read in some of the signature bytes */
   if (fread(buf, 1, PNG_BYTES_TO_CHECK, fp) != PNG_BYTES_TO_CHECK)
     error=NOT_ENOUGH_SIG_BYTES;

  /* Compare the first PNG_BYTES_TO_CHECK bytes of the signature.*/
   if (png_sig_cmp(buf, (png_size_t)0, PNG_BYTES_TO_CHECK))
     error=NOT_A_PNG_FILE;
  
   if (error)
     {
       return(error);
     }
   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also supply the
    * the compiler header file version, so that we know if the application
    * was compiled with a compatible version of the library.  REQUIRED
    */
   png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

   if (png_ptr == NULL)
   {
      fclose(fp);
      error=PROBLEM_READ_STRUCT;
      return (error);
   }

   /* Allocate/initialize the memory for image information.  REQUIRED. */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL)
   {
      fclose(fp);
      png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
      error=PROBLEM_INFO_STRUCT;
      return (error);
   }

   /* Set error handling if you are using the setjmp/longjmp method (this is
    * the normal method of doing things with libpng).  REQUIRED unless you
    * set up your own error handlers in the png_create_read_struct() earlier.
    */

   if (setjmp(png_jmpbuf(png_ptr)))
   {
      /* Free all of the memory associated with the png_ptr and info_ptr */
      png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
      fclose(fp);
      /* If we get here, we had a problem reading the file */
      error=PROBLEM_JMPBUF;
      return (error);
   }

   /* Set up the input control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* If we have already read some of the signature */
   png_set_sig_bytes(png_ptr, sig_read);

   /*
    * If you have enough memory to read in the entire image at once,
    * and you need to specify only transforms that can be controlled
    * with one of the PNG_TRANSFORM_* bits (this presently excludes
    * dithering, filling, setting background, and doing gamma
    * adjustment), then you can read the entire image (including
    * pixels) into the info structure with this call:
    */
   png_transforms=PNG_TRANSFORM_STRIP_16;
   png_read_png(png_ptr, info_ptr, png_transforms, NULL);
   /* At this point you have read the entire image */

   png_get_IHDR(png_ptr, info_ptr, &width, &height,
		&bit_depth, &color_type, &interlace_type,
		&compression_type, &filter_type);

   //printf("width %d\n",(int)width);
   //printf("height %d\n",(int)height);
   //printf("bit_depth %d\n",bit_depth);
   //printf("color_type %d\n",color_type);
   //printf("interlace_type %d\n",interlace_type);
   //printf("compression_type %d\n",compression_type);
   //printf("filter_type %d\n",filter_type);
   
   row_bytes=png_get_rowbytes(png_ptr,info_ptr);

   //printf("Number of bytes required to hold a row %d\n",row_bytes);

   //row_pointers = png_malloc(png_ptr, height*sizeof(png_bytep));

   row_pointers = png_get_rows(png_ptr, info_ptr);

   //   for (i=0;i<height;i++) printf("%d %d %d \n",row_pointers[i][0],row_pointers[i][100],row_pointers[i][200]);
   
   fclose(fp);
   return(IMAGE_OK);
}

// returns the number of nodes = number of pixels -1
guint32 png_image::get_nx()
{
  return((guint32) width);
}
// returns the number of nodes = number of pixels -1
guint32 png_image::get_ny()
{
  return((guint32) height);
}
guint32 png_image::get_error()
{
  return(error);
}

guint32 png_image::get_pixel(guint32 x, guint32 y)
{

  guint32 pixel;

  pixel=row_pointers[y][x*3];
  pixel|=(row_pointers[y][x*3+1] <<8);
  pixel|=(row_pointers[y][x*3+2] <<16);

  return(pixel);
}


#define OK 0
#define ERROR 1
/* write a png file */
guint32 png_image::write_image(char *file_name, guint32 x, guint32 y, png_bytep image, guint32 bytes_per_pixel)
{
   png_uint_32 k;//, height, width;
   //int bytes_per_pixel=3;
   //png_bytep image;
   int bit_depth;

   height=y;
   width=x;

   row_pointers=new png_bytep[height];
   for (k = 0; k < height; k++)
     row_pointers[k] = image + k*width*bytes_per_pixel;

   bit_depth=8;


   /* open the file */
   fp = fopen(file_name, "wb");
   if (fp == NULL)
      return (ERROR);

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also check that
    * the library version is compatible with the one used at compile time,
    * in case we are using dynamically linked libraries.  REQUIRED.
    */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL, NULL, NULL);

   if (png_ptr == NULL)
   {
      fclose(fp);
      return (ERROR);
   }

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL)
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
      return (ERROR);
   }
   //printf("done");

   //printf("error handling\n");

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr)))
   {
      /* If we get here, we had a problem reading the file */
      fclose(fp);
      png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
      return (ERROR);
   }
   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);
   //printf("done\n");

   /* This is the hard way */

   /* Set the image information here.  Width and height are up to 2^31,
    * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
    * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
    * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
    * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
    * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
    * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
    */
   //printf("IHDR\n");

   png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, PNG_COLOR_TYPE_RGB,
      PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

   /* Write the file header information.  REQUIRED */

   png_write_info(png_ptr, info_ptr);

   png_write_image(png_ptr, row_pointers);

   png_write_end(png_ptr, info_ptr);

   delete[] row_pointers;
   delete[] image;

   png_destroy_write_struct(&png_ptr, &info_ptr);

   /* close the file */
   fclose(fp);

   /* that's it */
   return (OK);
}

