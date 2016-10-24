
/* example.c - an example of using libpng */

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

#include "png.h"

 /* The png_jmpbuf() macro, used in error handling, became available in
  * libpng version 1.0.6.  If you want to be able to run your code with older
  * versions of libpng, you must define the macro yourself (but only if it
  * is not already defined by libpng!).
  */

#ifndef png_jmpbuf
#  define png_jmpbuf(png_ptr) ((png_ptr)->jmpbuf)
#endif


#define OK 0
#define ERROR 1
/* write a png file */
int write_png(char *file_name) /*, height, width); */
{
   FILE *fp;
   png_structp png_ptr;
   png_infop info_ptr;
   png_colorp palette;
   png_uint_32 k, height, width;
   int bytes_per_pixel;
/*    png_bytep row_pointers[height]; */
/*     png_byte image[height][width*bytes_per_pixel]; */
   png_bytep row_pointers[512];
   png_byte image[512][512*2];
   int true_red_bit_depth, true_green_bit_depth,true_blue_bit_depth;
   int bit_depth;
   int i,j;

   height=512;
   width=512;

   for (i=0;i<height;i++)
     {
     for (j=0;j<width*2;j++)
       {
	 image[i][j]=100;
       }
     }

   true_red_bit_depth=6;
   true_green_bit_depth=5;
   true_blue_bit_depth=5;

   bytes_per_pixel=2;
   height=512;
   width=512;

   bit_depth=16;

   printf("assignments\n");

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
/*     png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, */
/*        png_voidp user_error_ptr, user_error_fn, user_warning_fn); */
   printf("writestruct");

   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL);
   printf(" done\n");

   if (png_ptr == NULL)
   {
      fclose(fp);
      return (ERROR);
   }
   printf("infostruct\n");

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL)
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
      return (ERROR);
   }
   printf("done");

   printf("error handling\n");

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr)))
   {
      /* If we get here, we had a problem reading the file */
      fclose(fp);
      png_destroy_write_struct(&png_ptr,  (png_infopp)NULL);
/*        png_destroy_write_struct(&png_ptr, &info_ptr); */
      return (ERROR);
   }
   printf("done\n");

   printf("io\n");
   /* One of the following I/O initialization functions is REQUIRED */
   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);
   printf("done\n");

   /* This is the hard way */

   /* Set the image information here.  Width and height are up to 2^31,
    * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
    * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
    * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
    * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
    * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
    * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
    */
   printf("IHDR\n");

   png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, PNG_COLOR_TYPE_RGB,
      PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
   printf("done\n");

   /* set the palette if there is one.  REQUIRED for indexed-color images */
/*     palette = (png_colorp)png_malloc(png_ptr, 256 * sizeof (png_color)); */
   /* ... set palette colors ... */
/*     png_set_PLTE(png_ptr, info_ptr, palette, 256); */
   /* You must not free palette here, because png_set_PLTE only makes a link to
      the palette that you malloced.  Wait until you are about to destroy
      the png structure. */

   /* optional significant bit chunk */
   /* if we are dealing with a grayscale image then */
/*     sig_bit.gray = true_bit_depth; */
   /* otherwise, if we are dealing with a color image then */
/*     sig_bit.red = true_red_bit_depth; */
/*     sig_bit.green = true_green_bit_depth; */
/*     sig_bit.blue = true_blue_bit_depth; */
   /* if the image has an alpha channel then */
/*     sig_bit.alpha = true_alpha_bit_depth; */
/*     png_set_sBIT(png_ptr, info_ptr, sig_bit); */


   /* Optional gamma chunk is strongly suggested if you have any guess
    * as to the correct gamma of the image.
    */
/*     png_set_gAMA(png_ptr, info_ptr, gamma); */

   /* Optionally write comments into the image */
/*     text_ptr[0].key = "Title"; */
/*     text_ptr[0].text = "Mona Lisa"; */
/*     text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE; */
/*     text_ptr[1].key = "Author"; */
/*     text_ptr[1].text = "Leonardo DaVinci"; */
/*     text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE; */
/*     text_ptr[2].key = "Description"; */
/*     text_ptr[2].text = "<long text>"; */
/*     text_ptr[2].compression = PNG_TEXT_COMPRESSION_zTXt; */
#ifdef PNG_iTXt_SUPPORTED
/*     text_ptr[0].lang = NULL; */
/*     text_ptr[1].lang = NULL; */
/*     text_ptr[2].lang = NULL; */
#endif
/*     png_set_text(png_ptr, info_ptr, text_ptr, 3); */

   /* other optional chunks like cHRM, bKGD, tRNS, tIME, oFFs, pHYs, */
   /* note that if sRGB is present the gAMA and cHRM chunks must be ignored
    * on read and must be written in accordance with the sRGB profile */

   /* Write the file header information.  REQUIRED */
   printf("Info\n");

   png_write_info(png_ptr, info_ptr);
   printf("done\n");

   /* If you want, you can write the info in two steps, in case you need to
    * write your private chunk ahead of PLTE:
    *
    *   png_write_info_before_PLTE(write_ptr, write_info_ptr);
    *   write_my_chunk();
    *   png_write_info(png_ptr, info_ptr);
    *
    * However, given the level of known- and unknown-chunk support in 1.1.0
    * and up, this should no longer be necessary.
    */

   /* Once we write out the header, the compression type on the text
    * chunks gets changed to PNG_TEXT_COMPRESSION_NONE_WR or
    * PNG_TEXT_COMPRESSION_zTXt_WR, so it doesn't get written out again
    * at the end.
    */

   /* set up the transformations you want.  Note that these are
    * all optional.  Only call them if you want them.
    */

   /* invert monochrome pixels */
   /*     png_set_invert_mono(png_ptr); */

   /* Shift the pixels up to a legal bit depth and fill in
    * as appropriate to correctly scale the image.
    */
/*     png_set_shift(png_ptr, &sig_bit); */

   /* pack pixels into bytes */
/*     png_set_packing(png_ptr); */

   /* swap location of alpha bytes from ARGB to RGBA */
/*     png_set_swap_alpha(png_ptr); */

   /* Get rid of filler (OR ALPHA) bytes, pack XRGB/RGBX/ARGB/RGBA into
    * RGB (4 channels -> 3 channels). The second parameter is not used.
    */
/*     png_set_filler(png_ptr, 0, PNG_FILLER_BEFORE); */

   /* flip BGR pixels to RGB */
/*     png_set_bgr(png_ptr); */

   /* swap bytes of 16-bit files to most significant byte first */
/*     png_set_swap(png_ptr); */

   /* swap bits of 1, 2, 4 bit packed pixel formats */
/*     png_set_packswap(png_ptr); */

   /* turn on interlace handling if you are not using png_write_image() */
/*     if (interlacing) */
/*        number_passes = png_set_interlace_handling(png_ptr); */
/*     else */
/*        number_passes = 1; */

   /* The easiest way to write the image (you may have a different memory
    * layout, however, so choose what fits your needs best).  You need to
    * use the first method if you aren't handling interlacing yourself.
    */
   printf("pointers\n");

   for (k = 0; k < height; k++)
     row_pointers[k] = &image[0][0] + k*width*bytes_per_pixel;
   printf("done\n");


   /* One of the following output methods is REQUIRED */
   /* write out the entire image data in one call */
   printf("write image\n");

   png_write_image(png_ptr, row_pointers);

   printf("done\n");

   /* You can write optional chunks like tEXt, zTXt, and tIME at the end
    * as well.  Shouldn't be necessary in 1.1.0 and up as all the public
    * chunks are supported and you can use png_set_unknown_chunks() to
    * register unknown chunks into the info structure to be written out.
    */

   /* It is REQUIRED to call this to finish writing the rest of the file */
   printf("write end\n");

   png_write_end(png_ptr, info_ptr);

   printf("done\n");

   /* If you png_malloced a palette, free it here (don't free info_ptr->palette,
      as recommended in versions 1.0.5m and earlier of this example; if
      libpng mallocs info_ptr->palette, libpng will free it).  If you
      allocated it with malloc() instead of png_malloc(), use free() instead
      of png_free(). */
/*     png_free(png_ptr, palette); */
/*     palette=NULL; */

   /* Similarly, if you png_malloced any data that you passed in with
      png_set_something(), such as a hist or trans array, free it here,
      when you can be sure that libpng is through with it. */
 /*    png_free(png_ptr, trans); */
/*     trans=NULL; */

   printf("destroy\n");
   /* clean up after the write, and free any memory allocated */
   png_destroy_write_struct(&png_ptr, &info_ptr);
   printf("done\n");

   /* close the file */
   fclose(fp);

   /* that's it */
   return (OK);
}

int main (int argc, char *argv[])
{
  write_png("Testimage.png");
  return(0);
}
