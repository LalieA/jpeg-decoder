#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "../include/huffman.h"
#include "../include/jpeg_reader.h"

/// Structure containing the information contained in the header of the file
struct jpeg_desc {
  struct bitstream *stream;
  char *filename;

  uint16_t img_width;
  uint16_t img_height;

  /* Encoding infos */
  uint8_t is_progressive;
  uint8_t successive_approx;

  /* Spectral selection indexes */
  uint8_t Ss;
  uint8_t Se;

  /* Sampling factors */
  uint8_t YV_factor;
  uint8_t YH_factor;
  uint8_t CbV_factor;
  uint8_t CbH_factor;
  uint8_t CrV_factor;
  uint8_t CrH_factor;

  /* Components & associated quant' tables indexes */
  uint8_t nb_components;
  uint8_t index_Y_component;
  uint8_t index_Cb_component;
  uint8_t index_Cr_component;
  uint8_t index_Y_component_quant;
  uint8_t index_Cb_component_quant;
  uint8_t index_Cr_component_quant;

  /* Scan indexes */
  uint8_t scans[3];

  uint8_t nb_quantization_tables;
  uint8_t **quantization_tables;

  /* Huffman AC tables */
  uint8_t nb_huffman_tables_ac;
  struct huff_table **huff_table_ac;

  uint8_t index_huff_Y_ac;
  uint8_t index_huff_Cb_ac;
  uint8_t index_huff_Cr_ac;

  /* Huffman DC tables */
  uint8_t nb_huffman_tables_dc;
  struct huff_table **huff_table_dc;

  uint8_t index_huff_Y_dc;
  uint8_t index_huff_Cb_dc;
  uint8_t index_huff_Cr_dc;
};

void bad_marker(int expected, int got) {
  fprintf(stderr, "ERROR: Invalid marker, expected %x, got %x\n", expected,
          got);
  abort();
}

void not_jfif() {
  fprintf(stderr, "ERROR: Not JFIF 1.1 format\n");
  abort();
}

void goto_next_marker(struct jpeg_desc *jpeg, uint32_t *bytes) {
  while (*bytes != 0xff && !bitstream_is_empty(jpeg->stream)) {
    bitstream_read(jpeg->stream, 8, bytes, false);
  }
  bitstream_read(jpeg->stream, 8, bytes, false);
}

void read_sof(struct jpeg_desc *jpeg, uint32_t *bytes) {
  bitstream_read(jpeg->stream, 16, bytes, false); // skip section len
  bitstream_read(jpeg->stream, 8, bytes,
                 false); // check bit precision (always 8 for baseline)
  if (*bytes != 8) {
    fprintf(stderr, "ERROR: Bad bit precision, expected %x, got %x\n", 8,
            *bytes);
    abort();
  }

  bitstream_read(jpeg->stream, 16, bytes, false);
  jpeg->img_height = *bytes;
  bitstream_read(jpeg->stream, 16, bytes, false);
  jpeg->img_width = *bytes;

  bitstream_read(jpeg->stream, 8, bytes, false);
  jpeg->nb_components = *bytes;

  bitstream_read(jpeg->stream, 8, bytes, false);
  jpeg->index_Y_component = *bytes;
  bitstream_read(jpeg->stream, 4, bytes, false);
  jpeg->YH_factor = *bytes;
  bitstream_read(jpeg->stream, 4, bytes, false);
  jpeg->YV_factor = *bytes;
  bitstream_read(jpeg->stream, 8, bytes, false);
  jpeg->index_Y_component_quant = *bytes;

  if (jpeg->nb_components != 1) {
    bitstream_read(jpeg->stream, 8, bytes, false);
    jpeg->index_Cb_component = *bytes;
    bitstream_read(jpeg->stream, 4, bytes, false);
    jpeg->CbH_factor = *bytes;
    bitstream_read(jpeg->stream, 4, bytes, false);
    jpeg->CbV_factor = *bytes;
    bitstream_read(jpeg->stream, 8, bytes, false);
    jpeg->index_Cb_component_quant = *bytes;

    bitstream_read(jpeg->stream, 8, bytes, false);
    jpeg->index_Cr_component = *bytes;
    bitstream_read(jpeg->stream, 4, bytes, false);
    jpeg->CrH_factor = *bytes;
    bitstream_read(jpeg->stream, 4, bytes, false);
    jpeg->CrV_factor = *bytes;
    bitstream_read(jpeg->stream, 8, bytes, false);
    jpeg->index_Cr_component_quant = *bytes;
  }
  goto_next_marker(jpeg, bytes);
}

void read_app0(struct jpeg_desc *jpeg, uint32_t *bytes) {
  uint16_t len;

  bitstream_read(jpeg->stream, 16, bytes, false);
  if (*bytes != 0xffe0) {
    bad_marker(0xffe0, *bytes);
  };
  bitstream_read(jpeg->stream, 16, bytes, false);
  len = *bytes;
  bitstream_read(jpeg->stream, 8, bytes, false);
  if (*bytes != 'J') {
    not_jfif();
  }
  bitstream_read(jpeg->stream, 8, bytes, false);
  if (*bytes != 'F') {
    not_jfif();
  }
  bitstream_read(jpeg->stream, 8, bytes, false);
  if (*bytes != 'I') {
    not_jfif();
  }
  bitstream_read(jpeg->stream, 8, bytes, false);
  if (*bytes != 'F') {
    not_jfif();
  }
  bitstream_read(jpeg->stream, 8, bytes, false);
  if (*bytes != '\0') {
    not_jfif();
  }
  bitstream_read(jpeg->stream, 8, bytes, false);
  if (*bytes != 0x01) {
    not_jfif();
  }
  bitstream_read(jpeg->stream, 8, bytes, false);
  if (*bytes != 0x01) {
    not_jfif();
  }
  for (uint8_t i = 11; i < len; i++) {
    bitstream_read(jpeg->stream, 8, bytes, false);
  }
}

void read_dht(struct jpeg_desc *jpeg, uint32_t *bytes, uint16_t *nb_byte_read) {
  uint16_t len, cpt, index;

  bitstream_read(jpeg->stream, 16, bytes, false);
  len = *bytes - 2;
  for (cpt = 0; cpt < len; cpt++) {
    bitstream_read(jpeg->stream, 3, bytes, false);
    if (*bytes != 0) {
      fprintf(stderr, "ERROR: must be %x, got %x\n", 0, *bytes);
      abort();
    }

    bitstream_read(jpeg->stream, 1, bytes, false);
    enum sample_type type = *bytes;
    nb_byte_read = &index;

    if (type == DC) {
      jpeg->nb_huffman_tables_dc++;
      bitstream_read(jpeg->stream, 4, bytes, false);
      if (*bytes > 3) {
        fprintf(stderr, "ERROR: Must in 0..3, got %x\n", *bytes);
        abort();
      }
      index = *bytes;
      jpeg->huff_table_dc[index] =
          huffman_load_table(jpeg->stream, nb_byte_read);
    } else {
      jpeg->nb_huffman_tables_ac++;
      bitstream_read(jpeg->stream, 4, bytes, false);
      index = *bytes;
      jpeg->huff_table_ac[index] =
          huffman_load_table(jpeg->stream, nb_byte_read);
    }
    cpt += *nb_byte_read;
    goto_next_marker(jpeg, bytes);
  }
}

void read_sos(struct jpeg_desc *jpeg, uint32_t *bytes) {
  uint16_t len;
  bitstream_read(jpeg->stream, 16, bytes, false);
  len = *bytes;
  bitstream_read(jpeg->stream, 8, bytes, false);
  jpeg->nb_components = *bytes;
  if (len != 2 * *bytes + 6) {
    fprintf(stderr, "ERROR: Bad section length, expected %x, got %x\n",
            2 * *bytes + 6, len);
  }

  bitstream_read(jpeg->stream, 8, bytes, false);
  jpeg->scans[0] = *bytes;
  bitstream_read(jpeg->stream, 4, bytes, false);
  jpeg->index_huff_Y_dc = *bytes;
  bitstream_read(jpeg->stream, 4, bytes, false);
  jpeg->index_huff_Y_ac = *bytes;

  if (jpeg->nb_components != 1) {
    bitstream_read(jpeg->stream, 8, bytes, false);
    jpeg->scans[1] = *bytes;
    bitstream_read(jpeg->stream, 4, bytes, false);
    jpeg->index_huff_Cb_dc = *bytes;
    bitstream_read(jpeg->stream, 4, bytes, false);
    jpeg->index_huff_Cb_ac = *bytes;

    bitstream_read(jpeg->stream, 8, bytes, false);
    jpeg->scans[2] = *bytes;
    bitstream_read(jpeg->stream, 4, bytes, false);
    jpeg->index_huff_Cr_dc = *bytes;
    bitstream_read(jpeg->stream, 4, bytes, false);
    jpeg->index_huff_Cr_ac = *bytes;
  }
  bitstream_read(jpeg->stream, 8, bytes, false);
  if (!jpeg->is_progressive && *bytes != 0) {
    abort();
  }
  jpeg->Ss = *bytes;

  bitstream_read(jpeg->stream, 8, bytes, false);
  if (!jpeg->is_progressive && *bytes != 0x3f) {
    abort();
  }
  jpeg->Se = *bytes;

  bitstream_read(jpeg->stream, 8, bytes, false);
  if (!jpeg->is_progressive && *bytes != 0) {
    abort();
  }
  jpeg->successive_approx = *bytes;
  // goto_next_marker(jpeg, bytes);
}

void read_dqt(struct jpeg_desc *jpeg, uint32_t *bytes) {
  uint16_t len, index;

  bitstream_read(jpeg->stream, 16, bytes, false);
  len = *bytes - 2;
  jpeg->nb_quantization_tables += len / 65;

  for (uint8_t i = 0; i < len / 65; i++) {
    uint8_t *table = (uint8_t *)malloc(64 * sizeof(uint8_t));
    if (table == NULL) {
      abort();
    }
    bitstream_read(jpeg->stream, 4, bytes, false);
    bitstream_read(jpeg->stream, 4, bytes, false);
    index = *bytes;
    for (uint8_t j = 0; j < 64; j++) {
      bitstream_read(jpeg->stream, 8, bytes, false);
      table[j] = *bytes;
    }
    jpeg->quantization_tables[index] = table;
  }
  goto_next_marker(jpeg, bytes);
}

struct jpeg_desc *gen_jpeg_desc(const char *filename) {
  struct jpeg_desc *jpeg = (struct jpeg_desc *)malloc(sizeof(struct jpeg_desc));
  if (jpeg == NULL) {
    abort();
  }

  jpeg->stream = bitstream_create(filename);
  jpeg->filename = (char *)filename;

  jpeg->huff_table_ac =
      (struct huff_table **)malloc(4 * sizeof(struct huff_table *));
  if (jpeg->huff_table_ac == NULL) {
    abort();
  }
  jpeg->huff_table_dc =
      (struct huff_table **)malloc(4 * sizeof(struct huff_table *));
  if (jpeg->huff_table_dc == NULL) {
    abort();
  }
  jpeg->quantization_tables = (uint8_t **)malloc(4 * sizeof(uint8_t *));
  if (jpeg->quantization_tables == NULL) {
    abort();
  }

  jpeg->nb_quantization_tables = 0;
  jpeg->nb_huffman_tables_ac = 0;
  jpeg->nb_huffman_tables_dc = 0;
  jpeg->nb_components = 0;

  jpeg->CbV_factor = 0;
  jpeg->CbH_factor = 0;
  jpeg->CrV_factor = 0;
  jpeg->CrH_factor = 0;

  return jpeg;
}

/// Opens the JPEG file, reads the headers and store the information. Reading
/// stops after "SOS"
struct jpeg_desc *jpeg_first_read(const char *filename) {
  struct jpeg_desc *jpeg = gen_jpeg_desc(filename);

  uint16_t *nb_byte_read = NULL;
  uint32_t bytes;

  // SOI
  bitstream_read(jpeg->stream, 16, &bytes, false);
  if (bytes != 0xffd8) {
    bad_marker(0xffd8, bytes);
  }

  // APP0
  read_app0(jpeg, &bytes);

  // Go to first useful flag (skipping COM flags)
  while (bytes != 0xff && !bitstream_is_empty(jpeg->stream)) {
    bitstream_read(jpeg->stream, 8, &bytes, false);
  }
  bitstream_read(jpeg->stream, 8, &bytes, false);

  while (bytes != 0xda) {
    switch (bytes) {
    case 0xc0: // Baseline DCT (SOF0)
      jpeg->is_progressive = 0;
      read_sof(jpeg, &bytes);
      break;

    case 0xc2: // Progressive DCT (SOF2)
      jpeg->is_progressive = 1;
      read_sof(jpeg, &bytes);
      break;

    case 0xc4: // Define Huffman Tables (DHT)
      read_dht(jpeg, &bytes, nb_byte_read);
      break;

    case 0xdb: // Define Quantization tables (DQT)
      read_dqt(jpeg, &bytes);
      break;

    default:
      fprintf(stderr, "ERROR: Unknown flag %x\n", bytes);
      abort();
    }
  }
  read_sos(jpeg, &bytes); // Start of Scan (SOS)
  // bitstream_read(jpeg->stream, 8, &bytes, false);
  return jpeg;
}

void jpeg_read(struct jpeg_desc *jpeg) {
  uint16_t *nb_byte_read = NULL;
  uint32_t bytes;

  bitstream_read(jpeg->stream, 8, &bytes, false); // 0xff
  bitstream_read(jpeg->stream, 8, &bytes, false); // Flag
  while (bytes != 0xda) {
    switch (bytes) {
    case 0xc4: // Define Huffman Tables (DHT)
      read_dht(jpeg, &bytes, nb_byte_read);
      break;

    default:
      fprintf(stderr, "ERROR: Unknown flag %x\n", bytes);
      abort();
    }
  }
  read_sos(jpeg, &bytes); // Start of Scan (SOS)
}

/// Closes the file and frees all the required memory.
void jpeg_close(struct jpeg_desc *jpeg) {
  bitstream_close(jpeg->stream);

  for (uint8_t index = 0; index < jpeg_get_nb_quantization_tables(jpeg);
       index++) {
    free(jpeg->quantization_tables[index]);
  }
  for (uint8_t index = 0; index < jpeg_get_nb_huffman_tables(jpeg, AC);
       index++) {
    huffman_free_table(jpeg->huff_table_ac[index]);
  }
  for (uint8_t index = 0; index < jpeg_get_nb_huffman_tables(jpeg, DC);
       index++) {
    huffman_free_table(jpeg->huff_table_dc[index]);
  }

  free(jpeg->quantization_tables);
  free(jpeg->huff_table_ac);
  free(jpeg->huff_table_dc);
  free(jpeg);
}

/// Returns the filename of the open image.
char *jpeg_get_filename(const struct jpeg_desc *jpeg) { return jpeg->filename; }

/// Access to stream, placed just at the beginning of the scan raw data
struct bitstream *jpeg_get_bitstream(const struct jpeg_desc *jpeg) {
  return jpeg->stream;
}

// from DQT
/// Returns the number of quantization tables.
uint8_t jpeg_get_nb_quantization_tables(const struct jpeg_desc *jpeg) {
  return jpeg->nb_quantization_tables;
}

/// Returns the address of the i-th quantization table.
uint8_t *jpeg_get_quantization_table(const struct jpeg_desc *jpeg,
                                     uint8_t index) {
  return jpeg->quantization_tables[index];
}

// from DHT
/// Returns the number of Huffman tables.
uint8_t jpeg_get_nb_huffman_tables(const struct jpeg_desc *jpeg,
                                   enum sample_type acdc) {
  if (acdc == AC) {
    return jpeg->nb_huffman_tables_ac;
  }
  return jpeg->nb_huffman_tables_dc;
}

/// Returns the address of the i-th Huffman table of type "acdc".
struct huff_table *jpeg_get_huffman_table(const struct jpeg_desc *jpeg,
                                          enum sample_type acdc,
                                          uint8_t index) {
  if (acdc == AC) {
    return jpeg->huff_table_ac[index];
  }
  return jpeg->huff_table_dc[index];
}

// from Frame Header SOF0 / SOF2
/// Returns 1 if a progressive encoding was used.
uint8_t is_progressive(const struct jpeg_desc *jpeg) {
  return jpeg->is_progressive;
}

/// Returns the image size in the direction "dir".
uint16_t jpeg_get_image_size(const struct jpeg_desc *jpeg, enum direction dir) {
  if (dir == H) {
    return jpeg->img_width;
  }
  return jpeg->img_height;
}

/// Returns the number of component of the image.
uint8_t jpeg_get_nb_components(const struct jpeg_desc *jpeg) {
  return jpeg->nb_components;
}

/// Returns the identifier of of the i-th component.
uint8_t jpeg_get_frame_component_id(const struct jpeg_desc *jpeg,
                                    uint8_t frame_comp_index) {
  switch (frame_comp_index) {
  case 0:
    return jpeg->index_Y_component;
  case 1:
    return jpeg->index_Cb_component;
  case 2:
    return jpeg->index_Cr_component;
  }
  return 0;
}

/// Returns the sampling factor of of the i-th component in the direction "dir".
uint8_t jpeg_get_frame_component_sampling_factor(const struct jpeg_desc *jpeg,
                                                 enum direction dir,
                                                 uint8_t frame_comp_index) {
  switch (frame_comp_index) {
  case 0:
    if (dir == H) {
      return jpeg->YH_factor;
    } else {
      return jpeg->YV_factor;
    }
  case 1:
    if (dir == H) {
      return jpeg->CbH_factor;
    } else {
      return jpeg->CbV_factor;
    }
  case 2:
    if (dir == H) {
      return jpeg->CrH_factor;
    } else {
      return jpeg->CrV_factor;
    }
  }
  return 0;
}

/// Returns the index of of the i-th component.
uint8_t jpeg_get_frame_component_quant_index(const struct jpeg_desc *jpeg,
                                             uint8_t frame_comp_index) {
  switch (frame_comp_index) {
  case 0:
    return jpeg->index_Y_component_quant;
  case 1:
    return jpeg->index_Cb_component_quant;
  case 2:
    return jpeg->index_Cr_component_quant;
  }
  return 0;
}

// from Scan Header SOS
/// Returns the identifier of the i-th component read in the scan.
uint8_t jpeg_get_scan_component_id(const struct jpeg_desc *jpeg,
                                   uint8_t scan_comp_index) {
  return jpeg->scans[scan_comp_index];
}

/// Returns the index of the Huffman table (DC or AC) of the i-th component read
/// in the scan.
uint8_t jpeg_get_scan_component_huffman_index(const struct jpeg_desc *jpeg,
                                              enum sample_type acdc,
                                              uint8_t scan_comp_index) {
  switch (scan_comp_index) {
  case 0:
    if (acdc == AC) {
      return jpeg->index_huff_Y_ac;
    } else {
      return jpeg->index_huff_Y_dc;
    }
  case 1:
    if (acdc == AC) {
      return jpeg->index_huff_Cb_ac;
    } else {
      return jpeg->index_huff_Cb_dc;
    }
  case 2:
    if (acdc == AC) {
      return jpeg->index_huff_Cr_ac;
    } else {
      return jpeg->index_huff_Cr_dc;
    }
  }
  return 0;
}

uint8_t jpeg_get_Se(const struct jpeg_desc *jpeg) { return jpeg->Se; }

uint8_t jpeg_get_Ss(const struct jpeg_desc *jpeg) { return jpeg->Ss; }

uint8_t jpeg_get_successive_approx(const struct jpeg_desc *jpeg) {
  return jpeg->successive_approx;
}
