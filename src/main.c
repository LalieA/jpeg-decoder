#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "../include/jpeg_decode.h"
#include "../include/ppm_encode.h"

/// Call this function to give help.
void print_help(const char *exe_name) {
  printf(
      " \n"
      "JPEG to PGM/PPM image encoder\n"
      "   Decodes jpeg images and re-encodes them in portable pixmap format.\n"
      "   Make sure to pass only jpeg image files in arguments :)\n"
      "   Options:\n"
      "        -h, --help: shows this help\n"
      "        -v, --verbose: enables verbose\n"
      "        -o, --output: sets custom output for decompressed image\n"
      "   Usage:\n"
      "       %s -h | --help\n"
      "       %s -v | --verbose\n"
      "       %s filename_1.jpeg\n"
      "       %s filename_1.jpeg filename_2.jpg\n"
      "       %s filename_1.jpeg [-o | --output] output_filename.ppm\n"
      " \n"
      "Made with love by random Ensimag students <3\n"
      " \n",
      exe_name, exe_name, exe_name, exe_name, exe_name);
}

/// Call this function to print the correct usage.
void print_usage(const char *exe_name) {
  fprintf(stderr, "Usage: %s [-h] [-v] <file1> <file2> [-o <output2>]...\n",
          exe_name);
  printf("Re-run with \"-h\" to get help !\n");
}

int main(int argc, char **argv) {
  if (argc < 2) {
    /* Si y'a pas au moins un argument en ligne de commandes, on
     * boude. */
    print_usage(argv[0]);
    return EXIT_FAILURE;
  }

  char *custom_name;
  char *filename;
  u_int8_t verbose = 0;

  for (uint8_t i = 1; i < argc; i++) {
    if ((strcmp(argv[i], "-v") == 0) || strcmp(argv[i], "--verbose") == 0) {
      verbose = 1;
    }
  }

  for (uint8_t i = 1; i < argc; i++) {
    custom_name = NULL;
    if ((strcmp(argv[i], "-h") == 0) || strcmp(argv[i], "--help") == 0) {
      print_help(argv[0]);
      return EXIT_FAILURE;
    }
    if ((strcmp(argv[i], "-v") == 0) || strcmp(argv[i], "--verbose") == 0) {
      if (i < argc - 1) {
        i++;
      } else {
        return EXIT_SUCCESS;
      }
    } else {
      if (!is_jpeg(argv[i])) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
      };
    }

    /* On récupère le nom du fichier courant */
    filename = argv[i];

    if ((i < argc - 1) && (strcmp(argv[i + 1], "-o") == 0 ||
                           strcmp(argv[i + 1], "--output") == 0)) {
      i++;
      i++;
      custom_name = argv[i];
    }

    /* On créé une structure de fichier jpeg */
    struct jpeg_desc *jdesc = jpeg_first_read(filename);

    /* On récupère le flux des données brutes a partir du descripteur jpeg */
    struct bitstream *stream = jpeg_get_bitstream(jdesc);

    /* On décode le fichier courant et le réencode au format pgm/ppm */
    decode_jpeg(jdesc, stream, filename, init_new_file, jpeg_to_ppm,
                custom_name, verbose);

    /* Nettoyage de printemps : close_jpeg ferme aussi le bitstream
     * (voir Annexe C du sujet). */
    jpeg_close(jdesc);
  }

  /* On se congratule. */
  return EXIT_SUCCESS;
}
