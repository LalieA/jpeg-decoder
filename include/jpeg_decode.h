#ifndef __JPEG_DECODE_H__
#define __JPEG_DECODE_H__

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "bitstream.h"
#include "jpeg_reader.h"

struct block_info;
struct jpeg_work;

/// Structure that will contain a line of MCU stored in RGB.
struct mcu_line {
  uint8_t **mcu_array; // If gray scale, then a MCU is just one "layer". But if
                       // RGB, a mcu is represented by 3 juxtaposed layers
  uint8_t mcu_height;
  uint8_t mcu_width;
  uint16_t cur_pos; // index of the next MCU to consider
  uint16_t array_width;
};

/// Returns 1 of filename corresponds to a jpeg/jpg filename, 0 otherwise.
extern uint8_t is_jpeg(const char *filename);

/// Decodes a jpeg file and calls "write_func" to write RGB lines of MCU into
/// another format.
extern void decode_jpeg(const struct jpeg_desc *jdesc, struct bitstream *stream,
                        const char *filename,
                        char *(*init_func)(const struct jpeg_desc *jdesc,
                                           const char *filename, uint16_t width,
                                           uint16_t height, char *custom_name),
                        void (*write_func)(const struct jpeg_desc *jdesc,
                                           const struct mcu_line *mcu_line,
                                           const char *new_filename,
                                           uint16_t width, uint16_t height,
                                           uint16_t *y_written),
                        char *custom_name, uint8_t verbose);

#endif
