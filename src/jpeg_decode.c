#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/huffman.h"
#include "../include/jpeg_decode.h"
#include "../include/jpeg_reader.h"
#include "../include/verbose.h"

/// Generates and returns a block of size height*width with elements of type
/// "type_t".
#define gen_block(width, height, type_t)                                       \
  ({                                                                           \
    type_t **block = (type_t **)malloc(height * sizeof(type_t *));             \
    if (block == NULL) {                                                       \
      abort();                                                                 \
    }                                                                          \
    for (uint8_t i = 0; i < height; i++) {                                     \
      block[i] = (type_t *)malloc(width * sizeof(type_t));                     \
      if (block[i] == NULL) {                                                  \
        abort();                                                               \
      }                                                                        \
    }                                                                          \
    block;                                                                     \
  })

/// Frees the block given in argument.
#define free_block(block, height)                                              \
  ({                                                                           \
    for (uint8_t i = 0; i < height; i++) {                                     \
      free(*(block + i));                                                      \
    }                                                                          \
    free(block);                                                               \
  })

/// Some constants used to build an array.
#define A 0.4903926402016152
#define B 0.46193976625564337
#define C 0.4157348061512726
#define D 0.27778511650980114
#define E 0.19134171618254492
#define F 0.09754516100806417
#define G 0.35355339059327373

/// An array used in the IDCT step.
const float T[8][8] = {
    {G, G, G, G, G, G, G, G},     {A, C, D, F, -F, -D, -C, -A},
    {B, E, -E, -B, -B, -E, E, B}, {C, -F, -A, -D, D, A, F, -C},
    {G, -G, -G, G, G, -G, -G, G}, {D, -A, F, C, -C, -F, A, -D},
    {E, -B, B, -E, -E, B, -B, E}, {F, -D, C, -A, A, -C, D, -F}};

/// The pi value
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/// The sqrt(2) value
#ifndef SQRT_2
#define SQRT_2 1.414213562373095048801
#endif

#define BLOCK_VEC_LEN 64
#define BLOCK_LEN 8

/// This structure represents a YCbCr MCU. Each component is a uint8_t block.
struct ycbcr_mcu {
  uint8_t **blockY;
  uint8_t **blockCr;
  uint8_t **blockCb;

  uint8_t width;
  uint8_t height;
};

/// Structure used to store the block extracted from the bitstream and the
/// previous DC value for each component.
struct block_info {
  int16_t *block_vec; /// The block extracted represented by a vector
  int16_t dcY;        /// The previous DC value of the Y component
  int16_t dcCb;       /// The previous DC value of the Cb component
  int16_t dcCr;       /// The previous DC value of the Cr component
};

/// Used to store data that will be used and changed many times during the
/// process.
struct jpeg_work {
  float **Z_T;    /// 2D float array used in IDCT
  float **Tt_Z_T; /// 2D float array used in IDCT
  int16_t **block;
  uint8_t **spatial_block;
  struct ycbcr_mcu *ycbcr_mcu;
};

/// Structure containing the next color component to read and the number of the
/// same color component already read in the cycle.
struct component_info {
  enum color_component color_component;
  uint8_t index;
};

/// Returns 1 of filename corresponds to a jpeg/jpg filename, 0 otherwise.
uint8_t is_jpeg(const char *filename) {
  uint8_t res = 1;
  char *ext = strrchr(filename, '.');
  if (!ext || (strcmp(".jpeg", ext) != 0 && strcmp(".jpg", ext) != 0)) {
    res = 0;
  }
  return res;
}

/// Converts a number represented by its magnitude and index to a decimal number
int16_t mag_to_dec(int8_t mag, uint32_t index) {
  int16_t res;
  if (mag == 0) {
    return 0; // To avoid negative power
  }
  if ((index >> (mag - 1)) == 0) { // Then res must be negative
    res = -(1 << (mag)) + 1;
    res += (int16_t)index;
    return res;
  } else { // Then res is positive
    res = (int16_t)index - (1 << (mag - 1));
    res += (int32_t)1 << (mag - 1);
    return res;
  }
}

/// This function changes component_info so that for the next block read,
/// component_info->color_component contains the exact color component of the
/// block read.
void get_next_component(const struct jpeg_desc *jdesc,
                        struct component_info **component) {
  uint8_t nb_components = jpeg_get_nb_components(jdesc);

  if ((*component) == NULL) { // Initialisation
    *component = (struct component_info *)malloc(sizeof(struct component_info));
    if (*component == NULL) {
      abort();
    }
    (*component)->color_component = Y;
    (*component)->index = 0;
    return;
  }

  // If gray scale, then the next color component is always Y
  if (nb_components == 1) {
    (*component)->color_component = Y;
    (*component)->index = 0;
  }

  // If not gray scale, there are H_factor * V_factor blocks for each color of a
  // MCU
  else {
    uint8_t H_factor, V_factor; // Horizontal and Vertical sampling factors
    enum color_component old_component = (*component)->color_component;
    H_factor = jpeg_get_frame_component_sampling_factor(jdesc, H,
                                                        (uint8_t)old_component);
    V_factor = jpeg_get_frame_component_sampling_factor(jdesc, V,
                                                        (uint8_t)old_component);

    if ((*component)->index < H_factor * V_factor - 1) {
      // Then the next component is the same than the old one
      (*component)->index++;
    } else {
      uint8_t new_component = (uint8_t)old_component + 1;

      if (new_component == 3) {
        new_component =
            0; // The "fourth" component is Y (= restart at the begining)
      }
      (*component)->color_component = (enum color_component)new_component;
      (*component)->index = 0;
    }
  }
}

/// Returns the number of block in a line of MCU
uint32_t blocks_in_mcu_line(const struct jpeg_desc *jdesc) {
  uint16_t img_width = jpeg_get_image_size((struct jpeg_desc *)jdesc, H);

  uint8_t mcu_sf_H = jpeg_get_frame_component_sampling_factor(jdesc, H, 0);
  uint8_t mcu_sf_V = jpeg_get_frame_component_sampling_factor(jdesc, V, 0);

  uint32_t nb_mcu_line = img_width / (BLOCK_LEN * mcu_sf_H);
  if (img_width % (mcu_sf_H * BLOCK_LEN) != 0) {
    nb_mcu_line++;
  };

  uint8_t nb_block_mcu = mcu_sf_H * mcu_sf_V;

  mcu_sf_H = jpeg_get_frame_component_sampling_factor(jdesc, H, 1);
  mcu_sf_V = jpeg_get_frame_component_sampling_factor(jdesc, V, 1);
  nb_block_mcu = nb_block_mcu + mcu_sf_H * mcu_sf_V;

  mcu_sf_H = jpeg_get_frame_component_sampling_factor(jdesc, H, 2);
  mcu_sf_V = jpeg_get_frame_component_sampling_factor(jdesc, V, 2);
  nb_block_mcu = nb_block_mcu + mcu_sf_H * mcu_sf_V;

  return nb_mcu_line * nb_block_mcu;
}

/// Generates the structure "block_info" from the previous one. If the previous
/// is None, initialize it.
struct block_info *gen_block_info(struct block_info *previous_block_info) {
  struct block_info *block_info =
      (struct block_info *)malloc(sizeof(struct block_info));
  if (block_info == NULL) {
    abort();
  }
  block_info->block_vec = (int16_t *)malloc(sizeof(int16_t) * BLOCK_VEC_LEN);
  if (block_info->block_vec == NULL) {
    abort();
  }

  if (previous_block_info == NULL) {
    block_info->dcY = 0;
    block_info->dcCb = 0;
    block_info->dcCr = 0;
    return block_info;
  } else {
    block_info->dcY = previous_block_info->dcY;
    block_info->dcCb = previous_block_info->dcCb;
    block_info->dcCr = previous_block_info->dcCr;
    return block_info;
  }
}

/// Frees the structure "block_info" given in argument.
void free_block_info(struct block_info *block_info) {
  free(block_info->block_vec);
  free(block_info);
}

/// Reads a block from the bitstream and modify block_info->block_vec
/// Returns a uint8_t representing if the block read is the last of a line of
/// MCU in the image
uint8_t extract_block(const struct jpeg_desc *jdesc, struct bitstream *stream,
                      struct block_info *block_info,
                      const enum color_component next_component,
                      uint32_t *block_count, uint32_t blocks_in_mcu_line) {
  if (*block_count >= blocks_in_mcu_line) {
    return 1;
  }

  uint8_t Se, Ss, counter;
  Se = jpeg_get_Se(jdesc);
  Ss = jpeg_get_Ss(jdesc);

  (*block_count)++;

  // Get the Huffman tables
  struct huff_table *table_ac, *table_dc;
  if (next_component == Y) {
    table_ac = jpeg_get_huffman_table(jdesc, AC, 0);
    table_dc = jpeg_get_huffman_table(jdesc, DC, 0);
  } else {
    table_ac = jpeg_get_huffman_table(jdesc, AC, 1);
    table_dc = jpeg_get_huffman_table(jdesc, DC, 1);
  }

  // Get the DC value
  if (!is_progressive(jdesc) || (is_progressive(jdesc) && Se == 0)) {
    int8_t mag = huffman_next_value(table_dc, stream);

    uint32_t index;
    uint8_t nb_read = bitstream_read(stream, mag, &index, true);

    if ((uint32_t)nb_read != (uint8_t)mag) {
      printf("ERROR: can't read enough bits\n");
    }

    int16_t dc = mag_to_dec(mag, index);
    // printf("mag = %i, index = %i, dc = %i, ancien_dc = %i\n", mag, index, dc,
    // ancien_dc);

    int16_t real_value;
    if (next_component == Y) {
      real_value = dc + block_info->dcY;
      block_info->dcY = real_value;
    } else if (next_component == Cb) {
      real_value = dc + block_info->dcCb;
      block_info->dcCb = real_value;
    } else {
      real_value = dc + block_info->dcCr;
      block_info->dcCr = real_value;
    }
    block_info->block_vec[0] = real_value;
  }

  // Now let's read the 63 other values of the block
  counter = Se * is_progressive(jdesc) + 1;
  while ((is_progressive(jdesc) && counter <= Ss) ||
         (!is_progressive(jdesc) && counter < BLOCK_VEC_LEN)) {
    /* 4 possible cases
     0x00 : EndOfBlock
     0xF0 : 16 zeros
     0x?0 : error
     0xmn : m zeros and then a value != 0 with magnitude n
    */

    // Read 2 bytes in the bitstream
    int8_t rle_symbol = huffman_next_value(table_ac, stream);

    if (rle_symbol == 0) {
      while ((is_progressive(jdesc) && counter <= Ss) ||
             (!is_progressive(jdesc) && counter < BLOCK_VEC_LEN)) {
        block_info->block_vec[counter] = 0;
        counter++;
      }
      break;
    }

    else if ((uint8_t)rle_symbol == 0xF0) {
      for (uint8_t i = 0; i < 16; i++) {
        block_info->block_vec[counter + i] = 0;
      }
      counter += 16;
    }

    else if ((rle_symbol & 0x0F) == 0) {
      block_info->block_vec[counter] = 0;
      printf("Oups, something went wrong (invalid symbol) \n");
    }

    else {
      uint8_t zeros = (uint8_t)rle_symbol >> 4;
      for (int8_t i = 0; i < zeros; i++) {
        block_info->block_vec[counter + i] = 0;
      }
      counter += (uint8_t)zeros;

      uint8_t mag = (uint8_t)rle_symbol & 0x0F;

      // debug :
      if (mag < 1 || mag > 10) {
        printf("Oups, something went wrong (invalid mag) \n");
      }

      uint32_t index;
      uint8_t nb_read = bitstream_read(stream, mag, &index, true);
      if ((uint32_t)nb_read != (uint8_t)mag) {
        printf("Error : can't read enough bits\n");
      }
      block_info->block_vec[counter] = mag_to_dec(mag, index);
      counter++;
    }
  }
  return 0;
}

/// Returns a structure "ycbcr_mcu" with the given width and height. This
/// function is only used by "init_ycbcr".
struct ycbcr_mcu *gen_ycbcr_mcu(uint8_t width, uint8_t height) {
  struct ycbcr_mcu *mcu = (struct ycbcr_mcu *)malloc(sizeof(struct ycbcr_mcu));
  mcu->width = width;
  mcu->height = height;
  return mcu;
}

/// Initialise and returns a structure "ycbcr_mcu" with the correct width and
/// height (taken from jdesc).
struct ycbcr_mcu *init_ycbcr(const struct jpeg_desc *jdesc) {
  uint8_t YV = jpeg_get_frame_component_sampling_factor(jdesc, V, 0);
  uint8_t YH = jpeg_get_frame_component_sampling_factor(jdesc, H, 0);
  uint8_t h_scale = YH * BLOCK_LEN;
  uint8_t v_scale = YV * BLOCK_LEN;

  struct ycbcr_mcu *ycbcr_mcu = gen_ycbcr_mcu(h_scale, v_scale);
  ycbcr_mcu->blockY = gen_block(h_scale, v_scale, uint8_t);

  if (jpeg_get_nb_components(jdesc) == 3) {
    ycbcr_mcu->blockCb = gen_block(h_scale, v_scale, uint8_t);
    ycbcr_mcu->blockCr = gen_block(h_scale, v_scale, uint8_t);
  } else { // grey scale
    ycbcr_mcu->blockCb = NULL;
    ycbcr_mcu->blockCr = NULL;
  }
  return ycbcr_mcu;
};

/// Frees the mcu given in argument.
void free_ycbcr_mcu(struct ycbcr_mcu *mcu) {
  free_block(mcu->blockY, mcu->height);
  if (mcu->blockCr != NULL && mcu->blockCb != NULL) {
    free_block(mcu->blockCr, mcu->height);
    free_block(mcu->blockCb, mcu->height);
  }
  free(mcu);
}

/// Reorganizes a 1x64 vector-frequence domain (zig-zag encode) into a
/// BLOCK_LENxBLOCK_LEN block-frequence domain (inverse zig-zag encode).
void inverse_zig_zag(const int16_t *block_vec, int16_t **block) {
  uint16_t line = 0, col = 0;
  for (uint8_t i = 0; i < BLOCK_VEC_LEN; i++) {
    block[line][col] = block_vec[i];
    if ((line + col) % 2 == 0) { /* ASCENT */
      if (line == 0 && col + 1 != BLOCK_LEN) {
        col++;
      } // RIGHT SHIFT
      else if (col + 1 == BLOCK_LEN) {
        line++;
      } // DOWN SHIFT
      else {
        line--;
        col++;
      }
    } else { /* DESCENT */
      if (col == 0 && line + 1 != BLOCK_LEN) {
        line++;
      } // DOWN SHIFT
      else if (line + 1 == BLOCK_LEN) {
        col++;
      } // RIGHT SHIFT
      else {
        line++;
        col--;
      }
    }
  }
}

/// Inverse quantization of a 1x64 vector-frequence domain.
void inverse_quantization(const struct jpeg_desc *jdesc, int16_t *block_vec,
                          uint8_t index) {
  uint8_t *quant_table = jpeg_get_quantization_table(jdesc, index);
  for (uint8_t i = 0; i < BLOCK_VEC_LEN; i++) {
    *(block_vec + i) *= *(quant_table + i);
  }
}

/// Saturates a float between two values a and b.
void sat_float(float *float_to_sat, uint8_t a, uint8_t b) {
  if (*float_to_sat < a) {
    *float_to_sat = a;
  } else if (*float_to_sat > b) {
    *float_to_sat = b;
  }
}

/// Creates RGB values from ycbcr values.
void ycbcr_to_rgb(const struct ycbcr_mcu *ycbcr_mcu,
                  struct mcu_line *mcu_line) {
  uint8_t is_grey = ycbcr_mcu->blockCb == NULL && ycbcr_mcu->blockCr == NULL;
  uint16_t j_start = mcu_line->mcu_width * mcu_line->cur_pos;

  if (is_grey) {
    float gray;
    for (uint8_t i = 0; i < ycbcr_mcu->height; i++) {
      for (uint16_t j = 0; j < ycbcr_mcu->width; j++) {
        gray = ycbcr_mcu->blockY[i][j];
        sat_float(&gray, 0, 255);
        mcu_line->mcu_array[i][j + j_start] = (uint8_t)gray;
      }
    }
  }

  else {
    float red;
    float green, blue;
    uint16_t j_start_R = j_start * 3;
    uint16_t j_start_G = j_start_R + mcu_line->mcu_width;
    uint16_t j_start_B = j_start_R + mcu_line->mcu_width * 2;

    for (uint8_t i = 0; i < ycbcr_mcu->height; i++) {
      for (uint16_t j = 0; j < ycbcr_mcu->width; j++) {
        red =
            ycbcr_mcu->blockY[i][j] + 1.402 * (ycbcr_mcu->blockCr[i][j] - 128);
        green = ycbcr_mcu->blockY[i][j] -
                0.34414 * (ycbcr_mcu->blockCb[i][j] - 128) -
                0.71414 * (ycbcr_mcu->blockCr[i][j] - 128);
        blue =
            ycbcr_mcu->blockY[i][j] + 1.772 * (ycbcr_mcu->blockCb[i][j] - 128);
        sat_float(&red, 0, 255);
        sat_float(&green, 0, 255);
        sat_float(&blue, 0, 255);
        mcu_line->mcu_array[i][j + j_start_R] = (uint8_t)red;
        mcu_line->mcu_array[i][j + j_start_G] = (uint8_t)green;
        mcu_line->mcu_array[i][j + j_start_B] = (uint8_t)blue;
      }
    }
  }
}

/// Fast IDCT method.
void fast_idct(int16_t **block, uint8_t **spatial_block, float **Z_T,
               float **Tt_Z_T) {
  for (uint8_t i = 0; i < BLOCK_LEN; i++) {
    for (uint8_t j = 0; j < BLOCK_LEN; j++) {
      Z_T[i][j] = 0;
      for (uint8_t k = 0; k < BLOCK_LEN; k++) {
        Z_T[i][j] += ((float)block[i][k]) * T[k][j];
      }
    }
  }
  for (uint8_t i = 0; i < BLOCK_LEN; i++) {
    for (uint8_t j = 0; j < BLOCK_LEN; j++) {
      Tt_Z_T[i][j] = 128;
      for (uint8_t k = 0; k < BLOCK_LEN; k++) {
        Tt_Z_T[i][j] += T[k][i] * (float)Z_T[k][j];
      }
    }
  }
  for (uint8_t i = 0; i < BLOCK_LEN; i++) {
    for (uint8_t j = 0; j < BLOCK_LEN; j++) {
      sat_float(&Tt_Z_T[i][j], 0, 255);
      spatial_block[i][j] = (uint8_t)Tt_Z_T[i][j];
    }
  }
}

/// Returns a mcu with oversampled blocks depending on the chosen pattern.
void oversampling(const struct jpeg_desc *jpeg, struct ycbcr_mcu *mcu) {
  if (mcu->blockCb == NULL && mcu->blockCr == NULL) {
    return;
  }

  uint8_t h_sampling = jpeg_get_frame_component_sampling_factor(jpeg, H, 1);
  uint8_t v_sampling = jpeg_get_frame_component_sampling_factor(jpeg, V, 1);
  uint8_t h_sampling_y = jpeg_get_frame_component_sampling_factor(jpeg, H, 0);
  uint8_t v_sampling_y = jpeg_get_frame_component_sampling_factor(jpeg, V, 0);

  if (h_sampling == 1 && h_sampling_y == 2) {
    for (uint8_t line = 0; line < BLOCK_LEN; line++) {
      for (int8_t col = BLOCK_LEN - 1; col > -1; col--) {
        mcu->blockCb[line][col * 2] = mcu->blockCb[line][col];
        mcu->blockCb[line][col * 2 + 1] = mcu->blockCb[line][col];
        mcu->blockCr[line][col * 2] = mcu->blockCr[line][col];
        mcu->blockCr[line][col * 2 + 1] = mcu->blockCr[line][col];
      }
    }
  }

  if (v_sampling == 1 && v_sampling_y == 2) {
    for (uint8_t col = 0; col < BLOCK_LEN * h_sampling_y; col++) {
      for (int8_t line = BLOCK_LEN - 1; line > -1; line--) {
        mcu->blockCb[line * 2][col] = mcu->blockCb[line][col];
        mcu->blockCb[line * 2 + 1][col] = mcu->blockCb[line][col];
        mcu->blockCr[line * 2][col] = mcu->blockCr[line][col];
        mcu->blockCr[line * 2 + 1][col] = mcu->blockCr[line][col];
      }
    }
  }
}

/// Inserts a block in a mcu in its corresponding position (prepared for over
/// sampling).
void build_mcu(const struct jpeg_desc *jdesc,
               const enum color_component component, uint8_t mcu_index,
               struct ycbcr_mcu *ycbcr_mcu, const uint8_t **block) {
  uint8_t start_i, start_j;
  uint8_t YH = jpeg_get_frame_component_sampling_factor(jdesc, H, 0);
  if (component == Y) {
    start_i = BLOCK_LEN * (mcu_index / YH);
    start_j = BLOCK_LEN * (mcu_index % YH);
    // Copy the block in its new place
    for (uint8_t i = 0; i < BLOCK_LEN; i++) {
      for (uint8_t j = 0; j < BLOCK_LEN; j++) {
        ycbcr_mcu->blockY[i + start_i][j + start_j] = block[i][j];
      }
    }
    return;
  }

  uint8_t YV = jpeg_get_frame_component_sampling_factor(jdesc, V, 0);
  uint8_t CbH = jpeg_get_frame_component_sampling_factor(jdesc, H, 1);
  if (component == Cb) {
    // same than the Y component, but with "mcu_index - YV*YH"
    uint8_t new_mcu_index = mcu_index - YV * YH;
    start_i = BLOCK_LEN * (new_mcu_index / CbH);
    start_j = BLOCK_LEN * (new_mcu_index % CbH);

    for (uint8_t i = 0; i < BLOCK_LEN; i++) {
      for (uint8_t j = 0; j < BLOCK_LEN; j++) {
        ycbcr_mcu->blockCb[i + start_i][j + start_j] = block[i][j];
      }
    }
    return;
  }

  if (component == Cr) {
    uint8_t CbV = jpeg_get_frame_component_sampling_factor(jdesc, V, 1);
    uint8_t CrH = jpeg_get_frame_component_sampling_factor(jdesc, H, 2);
    uint8_t new_mcu_index = mcu_index - YV * YH - CbH * CbV;
    start_i = BLOCK_LEN * (new_mcu_index / CrH);
    start_j = BLOCK_LEN * (new_mcu_index % CrH);

    for (uint8_t i = 0; i < BLOCK_LEN; i++) {
      for (uint8_t j = 0; j < BLOCK_LEN; j++) {
        ycbcr_mcu->blockCr[i + start_i][j + start_j] = block[i][j];
      }
    }
  }
}

/// Returns the number of blocks in a MCU.
uint8_t get_blocks_in_mcu(const struct jpeg_desc *jdesc) {
  uint8_t YH = jpeg_get_frame_component_sampling_factor(jdesc, H, 0);
  uint8_t YV = jpeg_get_frame_component_sampling_factor(jdesc, V, 0);
  uint8_t CbH = jpeg_get_frame_component_sampling_factor(jdesc, H, 1);
  uint8_t CbV = jpeg_get_frame_component_sampling_factor(jdesc, V, 1);
  uint8_t CrH = jpeg_get_frame_component_sampling_factor(jdesc, H, 2);
  uint8_t CrV = jpeg_get_frame_component_sampling_factor(jdesc, V, 2);
  return (YH * YV) + (CbH * CbV) + (CrV * CrH);
}

/// Returns the number of lines of MCU in the image.
uint16_t get_nb_lines_of_mcu(const struct jpeg_desc *jdesc) {
  uint8_t YV = jpeg_get_frame_component_sampling_factor(jdesc, V, 0);
  uint16_t height = jpeg_get_image_size(jdesc, V);
  if (height % (YV * BLOCK_LEN) != 0) {
    return height / (8 * YV) + 1;
  }
  return height / (8 * YV);
}

/// Returns the number of column of MCU in the image.
uint16_t get_nb_column_of_mcu(const struct jpeg_desc *jdesc) {
  uint8_t YV = jpeg_get_frame_component_sampling_factor(jdesc, H, 0);
  uint16_t width = jpeg_get_image_size(jdesc, H);
  if (width % (YV * BLOCK_LEN) != 0) {
    return width / (8 * YV) + 1;
  }
  return width / (8 * YV);
}

/// Return an empty line of MCU.
struct mcu_line *gen_mcu_line(struct jpeg_desc *jdesc) {
  struct mcu_line *mcu_line = malloc(sizeof(struct mcu_line));
  mcu_line->mcu_width =
      8 * jpeg_get_frame_component_sampling_factor(jdesc, H, 0);
  mcu_line->mcu_height =
      8 * jpeg_get_frame_component_sampling_factor(jdesc, V, 0);

  uint16_t width = jpeg_get_nb_components(jdesc) * mcu_line->mcu_width *
                   get_nb_column_of_mcu(jdesc);

  uint8_t **line = malloc(mcu_line->mcu_height * sizeof(uint8_t *));
  if (line == NULL) {
    abort();
  }
  for (uint8_t i = 0; i < mcu_line->mcu_height; i++) {
    line[i] = (uint8_t *)malloc(width * sizeof(uint8_t));
    if (line[i] == NULL) {
      abort();
    }
  }
  mcu_line->mcu_array = line;
  mcu_line->array_width = width;
  return mcu_line;
}

/// Frees the given "mcu_line".
void free_mcu_line(struct mcu_line *mcu_line) {
  for (uint8_t i = 0; i < mcu_line->mcu_height; i++) {
    free(mcu_line->mcu_array[i]);
  }
  free(mcu_line->mcu_array);
  free(mcu_line);
}

/// Generates and returns the strucutre jpeg_work.
struct jpeg_work *initial_blocks(const struct jpeg_desc *jdesc) {
  float **Z_T = gen_block(BLOCK_LEN, BLOCK_LEN, float);
  float **Tt_Z_T = gen_block(BLOCK_LEN, BLOCK_LEN, float);
  int16_t **block = gen_block(BLOCK_LEN, BLOCK_LEN, int16_t);
  uint8_t **spatial_block = gen_block(BLOCK_LEN, BLOCK_LEN, uint8_t);
  struct ycbcr_mcu *ycbcr_mcu = init_ycbcr(jdesc);

  struct jpeg_work *jpeg_work = (struct jpeg_work *)malloc(sizeof(
      struct jpeg_work)); //{T, Z_T, Tt_Z_T, block, spatial_block, ycbcr_mcu};
  if (jpeg_work == NULL) {
    abort();
  }

  jpeg_work->block = block;
  jpeg_work->Z_T = Z_T;
  jpeg_work->Tt_Z_T = Tt_Z_T;
  jpeg_work->spatial_block = spatial_block;
  jpeg_work->ycbcr_mcu = ycbcr_mcu;

  return jpeg_work;
}

/// Frees the given "jpeg_work".
void free_jpeg_work(struct jpeg_work *jpeg_work) {
  free_block(jpeg_work->Z_T, BLOCK_LEN);
  free_block(jpeg_work->Tt_Z_T, BLOCK_LEN);
  free_block(jpeg_work->block, BLOCK_LEN);
  free_block(jpeg_work->spatial_block, BLOCK_LEN);
  free_ycbcr_mcu(jpeg_work->ycbcr_mcu);
  free(jpeg_work);
}

/// Decodes the given stream until it reads a full MCU line. The decoded data is
/// stored in the structure "mcu_line".
void jpeg_to_rgb_mcus(const struct jpeg_desc *jdesc, struct bitstream *stream,
                      struct block_info *block_info, struct jpeg_work *jwork,
                      struct mcu_line *mcu_line, uint8_t verbose) {
  struct component_info *next_component_index = NULL;

  uint32_t block_count = 0, block_in_mcu_line = blocks_in_mcu_line(jdesc);
  uint8_t index, block_index = 0, is_last_block = 0,
                 blocks_in_mcu = get_blocks_in_mcu(jdesc);

  while (!is_last_block) {
    get_next_component(jdesc, &next_component_index);
    is_last_block = extract_block(jdesc, stream, block_info,
                                  next_component_index->color_component,
                                  &block_count, block_in_mcu_line);

    if (!is_last_block) {
      if (verbose) {
        print_block_vec(block_info->block_vec)
      }

      // Index of the quantization tables
      if (next_component_index->color_component == Y) {
        index = 0;
      } else {
        index = 1;
      }

      inverse_quantization(jdesc, block_info->block_vec, index);
      if (verbose) {
        print_iquant(block_info->block_vec)
      }

      inverse_zig_zag(block_info->block_vec, jwork->block);
      if (verbose) {
        print_izz(jwork->spatial_block)
      }

      fast_idct(jwork->block, jwork->spatial_block, jwork->Z_T, jwork->Tt_Z_T);
      if (verbose) {
        print_idct(jwork->spatial_block)
      }

      build_mcu(jdesc, next_component_index->color_component, block_index,
                jwork->ycbcr_mcu, (const uint8_t **)jwork->spatial_block);

      block_index++;

      if (block_index == blocks_in_mcu) {
        oversampling(jdesc, jwork->ycbcr_mcu);
        if (verbose) {
          print_ycbcr_mcu(jwork->ycbcr_mcu,
                          next_component_index->color_component);
        }

        block_index = 0;
        ycbcr_to_rgb(jwork->ycbcr_mcu, mcu_line);
        mcu_line->cur_pos++;
      }
    }
  }
  free(next_component_index);
}

/// Decodes a jpeg file and calls "write_func" to write RGB lines of MCU into
/// another format.
void decode_jpeg(const struct jpeg_desc *jdesc, struct bitstream *stream,
                 const char *filename,
                 char *(*init_func)(const struct jpeg_desc *jdesc,
                                    const char *filename, uint16_t width,
                                    uint16_t height, char *custom_name),
                 void (*write_func)(const struct jpeg_desc *jdesc,
                                    const struct mcu_line *mcu_line,
                                    const char *new_filename, uint16_t width,
                                    uint16_t height, uint16_t *y_written),
                 char *custom_name, uint8_t verbose) {
  uint16_t width = jpeg_get_image_size(jdesc, H);
  uint16_t height = jpeg_get_image_size(jdesc, V);

  // First writes useful information to new file's header
  char *new_filename = init_func(jdesc, filename, width, height, custom_name);

  uint16_t nb_lines_of_mcu = get_nb_lines_of_mcu(jdesc);

  struct jpeg_work *jwork = initial_blocks(jdesc);
  struct block_info *block_info = gen_block_info(NULL);
  struct mcu_line *mcu_line = gen_mcu_line((struct jpeg_desc *)jdesc);

  uint16_t lines_mcu_read = 0; // To not read too much MCUs
  uint16_t y_written = 0;      // To handle y-axis cut
  uint8_t scans = 0;

  /* First scan */

  while (lines_mcu_read != nb_lines_of_mcu) {
    /* Decodes the stream to have a complete line of MCUs */
    mcu_line->cur_pos = 0;
    jpeg_to_rgb_mcus(jdesc, stream, block_info, jwork, mcu_line, verbose);

    /* Saves a line of MCUs to new file */
    write_func(jdesc, mcu_line, new_filename, width, height, &y_written);
    lines_mcu_read++;
  }
  scans++;

  /* Other scans */
  // while (jpeg_get_Se(jdesc) < BLOCK_VEC_LEN - 1 && scans < 10) {
  // scans++;
  // printf("Scan %d\n", scans);
  // jpeg_read((struct jpeg_desc *)jdesc);

  // while(lines_mcu_read != nb_lines_of_mcu) {
  //     /* Decodes the stream to have a complete line of MCUs */
  //     mcu_line->cur_pos = 0;
  //     jpeg_to_rgb_mcus(jdesc, stream, block_info, jwork, mcu_line);

  //     /* Saves a line of MCUs to new file */
  //     write_func(jdesc, mcu_line, new_filename, width, height, &y_written);
  //     lines_mcu_read++;
  // }
  // }

  free_block_info(block_info);
  free(new_filename);
  free_mcu_line(mcu_line);
  free_jpeg_work(jwork);
}
