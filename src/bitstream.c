#include "../include/bitstream.h"
#include <stdio.h>
#include <stdlib.h>

struct bitstream {
  FILE *file;
  bool is_empty;
  uint8_t buffer;
  uint8_t curr_bit;
};

/// Create a flux, starting at the begining of the file filename.
struct bitstream *bitstream_create(const char *filename) {
  FILE *f = fopen(filename, "rb");
  if (f == NULL) {
    fprintf(stderr, "ERROR: Cannot create bitstream for %s\n", filename);
    abort();
  }

  fseek(f, 0, SEEK_SET);

  struct bitstream *stream =
      (struct bitstream *)malloc(sizeof(struct bitstream));
  if (stream == NULL) {
    abort();
  }

  stream->file = f;
  stream->curr_bit = 0;

  if (fread(&stream->buffer, sizeof(uint8_t), 1, f) !=
      0) { // Bitstream not empty
    stream->is_empty = false;
    return stream;
  }
  // if empty, close bitstream return
  bitstream_close(stream);
  return NULL;
}

/// Close the file and free the steam.
void bitstream_close(struct bitstream *stream) {
  fclose(stream->file);
  free(stream);
}

void fill_buffer(struct bitstream *stream) {
  if (fread(&stream->buffer, sizeof(uint8_t), 1, stream->file) == 0) {
    stream->is_empty = true;
  }
  stream->curr_bit = 0;
}

/// Read "nb_bits" from the stream and write in do "dest".
uint8_t bitstream_read(struct bitstream *stream, uint8_t nb_bits,
                       uint32_t *dest, bool discard_byte_stuffing) {
  uint8_t nb_read_bits = 0, read_bit;
  *dest = 0;

  while (!stream->is_empty && nb_read_bits < nb_bits) {
    read_bit = stream->buffer >> (7 - stream->curr_bit) & 1;
    *dest = (*dest << 1) + read_bit;

    stream->curr_bit++;
    if (stream->curr_bit == 8) {
      if (stream->buffer == 0xff) { // Byte stuffing
        fill_buffer(stream);
        if (stream->buffer == 0x00 && discard_byte_stuffing) {
          fill_buffer(stream);
        }
      } else {
        fill_buffer(stream);
      }
    }
    nb_read_bits++;
  }
  return nb_read_bits;
}

/// Return True if the bitstream is empty, False otherwise.
bool bitstream_is_empty(struct bitstream *stream) { return stream->is_empty; }