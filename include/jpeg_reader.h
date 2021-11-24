#ifndef __JPEG_DESC_H__
#define __JPEG_DESC_H__

#include <stdbool.h>
#include <stdint.h>

#include "bitstream.h"

enum color_component { Y, Cb, Cr, NB_COLOR_COMPONENTS };

enum direction { H, V, NB_DIRECTIONS };

enum sample_type { DC, AC, NB_SAMPLE_TYPES };

struct jpeg_desc;

// #### general
/// Opens the JPEG file, reads the headers and store the information. Reading
/// stops after "SOS"
extern struct jpeg_desc *jpeg_first_read(const char *filename);

extern void jpeg_read(struct jpeg_desc *jpeg);

/// Closes the file and frees all the required memory.
extern void jpeg_close(struct jpeg_desc *jpeg);

/// Returns the filename of the open image.
extern char *jpeg_get_filename(const struct jpeg_desc *jpeg);

/// Access to stream, placed just at the beginning of the scan raw data
extern struct bitstream *jpeg_get_bitstream(const struct jpeg_desc *jpeg);

// ##### from DQT
/// Returns the number of quantization tables.
extern uint8_t jpeg_get_nb_quantization_tables(const struct jpeg_desc *jpeg);

/// Returns the address of the i-th quantization table.
extern uint8_t *jpeg_get_quantization_table(const struct jpeg_desc *jpeg,
                                            uint8_t index);

// ##### from DHT
/// Returns the number of Huffman tables.
extern uint8_t jpeg_get_nb_huffman_tables(const struct jpeg_desc *jpeg,
                                          enum sample_type acdc);
/// Returns the address of the i-th Huffman table of type "acdc".
extern struct huff_table *jpeg_get_huffman_table(const struct jpeg_desc *jpeg,
                                                 enum sample_type acdc,
                                                 uint8_t index);

// #### from Frame Header SOF0 / SOF2
/// Returns 1 if a progressive encoding was used.
extern uint8_t is_progressive(const struct jpeg_desc *jpeg);

/// Returns the image size in the direction "dir".
extern uint16_t jpeg_get_image_size(const struct jpeg_desc *jpeg,
                                    enum direction dir);

/// Returns the number of component of the image.
extern uint8_t jpeg_get_nb_components(const struct jpeg_desc *jpeg);

/// Returns the identifier of of the i-th component.
extern uint8_t jpeg_get_frame_component_id(const struct jpeg_desc *jpeg,
                                           uint8_t frame_comp_index);

/// Returns the sampling factor of of the i-th component in the direction "dir".
extern uint8_t jpeg_get_frame_component_sampling_factor(
    const struct jpeg_desc *jpeg, enum direction dir, uint8_t frame_comp_index);

/// Returns the index of of the i-th component.
extern uint8_t
jpeg_get_frame_component_quant_index(const struct jpeg_desc *jpeg,
                                     uint8_t frame_comp_index);

// #### from Scan Header SOS
/// Returns the identifier of the i-th component read in the scan.
extern uint8_t jpeg_get_scan_component_id(const struct jpeg_desc *jpeg,
                                          uint8_t scan_comp_index);

/// Returns the index of the Huffman table (DC or AC) of the i-th component read
/// in the scan.
extern uint8_t jpeg_get_scan_component_huffman_index(
    const struct jpeg_desc *jpeg, enum sample_type, uint8_t scan_comp_index);
extern uint8_t jpeg_get_Se(const struct jpeg_desc *jpeg);
extern uint8_t jpeg_get_Ss(const struct jpeg_desc *jpeg);
extern uint8_t jpeg_get_successive_approx(const struct jpeg_desc *jpeg);

#endif
