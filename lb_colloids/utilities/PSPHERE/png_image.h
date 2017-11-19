#ifndef PNG_IMAGE_000
#define PNG_IMAGE_000

#include <png.h>
#include <glib.h>

#define IMAGE_OK 0
#define FILE_NOT_FOUND 1
#define NOT_ENOUGH_SIG_BYTES 2
#define NOT_A_PNG_FILE 4
#define PROBLEM_READ_STRUCT 8
#define PROBLEM_INFO_STRUCT 16
#define PROBLEM_JMPBUF 32

class png_image
{
 protected:
  png_structp png_ptr;
  png_infop info_ptr;
  png_uint_32 width, height;
  FILE *fp;
  guint32 error;
  int bit_depth, color_type, interlace_type, compression_type, filter_type;

 public:
  png_bytepp row_pointers;
  png_image();
  ~png_image();
 public:
  guint32 read_image(char *filename);
  guint32 write_image(char *file_name, guint32 x, guint32 y, png_bytep image, guint32 bytes_per_pixel);

 public:
  guint32 get_nx();
  guint32 get_ny();
  guint32 get_error();

  guint32 get_pixel(guint32 x,guint32 y);
  // guint32 get_pixel_conf(guint32 x, guint32 y);

};


#endif
