#ifndef __PPM_ENCODE_H__
#define ___PPM_ENCODE_H__

#include "./jpeg_decode.h"
#include <stdint.h>

/// Write in the ppm file a line of MCU stored in the structure mcu_line.
extern void jpeg_to_ppm(const struct jpeg_desc *jdesc,
                        const struct mcu_line *mcu_line,
                        const char *new_filename, uint16_t width,
                        uint16_t height, uint16_t *y_written);

/// Creates a new file (replaces if already existing) and write the header. Also
/// returns the name of the new file.
extern char *init_new_file(const struct jpeg_desc *jdesc, const char *filename,
                           uint16_t width, uint16_t height, char *custom_name);

#endif
