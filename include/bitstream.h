#ifndef __BITSTREAM_H__
#define __BITSTREAM_H__

#include <stdbool.h>
#include <stdint.h>

struct bitstream;

/// Create a flux, starting at the begining of the file filename.
extern struct bitstream *bitstream_create(const char *filename);

/// Close the file and free the steam.
extern void bitstream_close(struct bitstream *stream);

/// Read "nb_bits" from the stream and write in do "dest".
extern uint8_t bitstream_read(struct bitstream *stream, uint8_t nb_bits,
                              uint32_t *dest, bool discard_byte_stuffing);

/// Return True if the bitstream is empty, False otherwise.
extern bool bitstream_is_empty(struct bitstream *stream);

#endif
