#ifndef __HUFFMAN_H__
#define __HUFFMAN_H__

#include "bitstream.h"
#include <stdbool.h>
#include <stdint.h>

struct huff_table;

/// Reads the stream and builds a Huffman table. At the beginning, the stream
/// must be 3 bytes after the flag DHT.
extern struct huff_table *huffman_load_table(struct bitstream *stream,
                                             uint16_t *nb_byte_read);

/// Returns the next value reached by using the Huffman table according to the
/// bits read in the stream.
extern int8_t huffman_next_value(struct huff_table *table,
                                 struct bitstream *stream);

/// Frees the memory used by the Huffman tables.
extern void huffman_free_table(struct huff_table *table);

#ifdef BLABLA
extern int8_t huffman_next_value_count(struct huff_table *table,
                                       struct bitstream *stream,
                                       uint8_t *nb_bits_read);
#endif
#endif
