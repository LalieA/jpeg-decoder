#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>

#include "../include/jpeg_decode.h"

/// Returns the filename of the new file the program is generating.
/// Builds it from the filename of the jpeg file and its extension.
char *get_new_filename(const char *filename, const char *extension,
                       char *custom_name) {
  if (custom_name != NULL) {
    filename = custom_name;
  }
  size_t f_len = strlen(filename);
  size_t e_len = strlen(extension);
  size_t new_len = f_len + e_len + 1;

  uint8_t old_e_len = 0;

  char curr = filename[f_len - 1];
  while (curr != '.' && old_e_len < f_len) {
    curr = filename[f_len - 1 - old_e_len];
    old_e_len++;
  }

  char *new_filename =
      (char *)malloc((f_len - old_e_len + new_len + 2) * sizeof(char));
  if (new_filename == NULL) {
    abort();
  }

  if (old_e_len == f_len) {
    old_e_len = 0;
  }
  strncpy(new_filename, filename, f_len - old_e_len);
  new_filename[f_len - old_e_len] = '\0';

  strcat(new_filename, ".");
  strcat(new_filename, extension);

  return new_filename;
}

/// Creates a new file (replaces if already existing) and write the header. Also
/// returns the name of the new file.
char *init_new_file(const struct jpeg_desc *jdesc, const char *filename,
                    uint16_t width, uint16_t height, char *custom_name) {
  char *extension;
  uint8_t is_grey = jpeg_get_nb_components(jdesc) != 3;

  if (is_grey) {
    extension = "pgm";
  } else {
    extension = "ppm";
  }

  char *new_filename = get_new_filename(filename, extension, custom_name);
  remove(new_filename);
  FILE *f = fopen(new_filename, "w");
  if (f == NULL) {
    abort();
  }

  if (is_grey) {
    fprintf(f, "P5\n");
  } else {
    fprintf(f, "P6\n");
  }

  fprintf(f, "%d %d\n", width, height);
  fprintf(f, "255\n");
  fclose(f);
  return new_filename;
}

/// Write in the ppm file a line of MCU stored in the structure mcu_line.
void jpeg_to_ppm(const struct jpeg_desc *jdesc, const struct mcu_line *mcu_line,
                 const char *new_filename, uint16_t width, uint16_t height,
                 uint16_t *y_written) {
  uint8_t is_grey = jpeg_get_nb_components(jdesc) != 3;

  FILE *f = fopen(new_filename, "a");

  if (is_grey) { // Then, we can just copy the mcus_array line by line.
    for (uint8_t i = 0; i < mcu_line->mcu_height && *y_written < height; i++) {
      (*y_written)++;
      for (uint16_t j = 0; j < width; j++) {
        fwrite(&(mcu_line->mcu_array[i][j]), sizeof(uint8_t), 1, f);
      }
    }
  }

  else {
    uint16_t x_written;
    uint16_t mcus_in_line = mcu_line->array_width / (3 * mcu_line->mcu_width);

    for (uint8_t line = 0; line < mcu_line->mcu_height && *y_written < height;
         line++) { // For all lines,
      (*y_written)++;
      x_written = 0;
      for (uint16_t mcu = 0; mcu < mcus_in_line; mcu++) { // we take each MCU
        for (uint8_t col = 0; col < mcu_line->mcu_width && x_written < width;
             col++) { // And write the values of
                      // the 3 colors for each column in the MCU
          x_written++;
          fwrite(
              &(mcu_line->mcu_array[line][col + mcu * 3 * mcu_line->mcu_width]),
              sizeof(uint8_t), 1, f);
          fwrite(
              &(mcu_line->mcu_array[line][col + mcu * 3 * mcu_line->mcu_width +
                                          mcu_line->mcu_width]),
              sizeof(uint8_t), 1, f);
          fwrite(
              &(mcu_line->mcu_array[line][col + mcu * 3 * mcu_line->mcu_width +
                                          2 * mcu_line->mcu_width]),
              sizeof(uint8_t), 1, f);
        }
      }
    }
  }
  fclose(f);
}
