#include "../include/jpeg_reader.h"
#include <stdio.h>

#define BLOCK_VEC_LEN 64
#define BLOCK_LEN 8

#define VERBOSE

#ifdef VERBOSE
/// Debug tool: print the component given
#define print_component(component)                                             \
  ({                                                                           \
    printf("\n** component ");                                                 \
    switch (component)                                                         \
    case Y:                                                                    \
      printf("Y");                                                             \
    break;                                                                     \
  case Cr:                                                                     \
    printf("Cr");                                                              \
    break;                                                                     \
  case Cb:                                                                     \
    printf("Cb");                                                              \
    break;                                                                     \
    printf("\n");                                                              \
  });
/// Debug tool: print the given block number
#define print_block_num(num) ({ printf("* bloc %d", num); });
/// Debug tool: print a block vector in hexa format after extration.
#define print_block_vec(block_vec)                                             \
  ({                                                                           \
    printf("\n[  bloc] ");                                                     \
    for (uint8_t i = 0; i < BLOCK_VEC_LEN; i++) {                              \
      printf("%x ", block_vec[i] & 0xffff);                                    \
    }                                                                          \
    printf("\n");                                                              \
  });
/// Debug tool: print a block vector in hexa format after iquant.
#define print_iquant(block_vec)                                                \
  ({                                                                           \
    printf("\n[iquant] ");                                                     \
    for (uint8_t i = 0; i < BLOCK_VEC_LEN; i++) {                              \
      printf("%x ", block_vec[i] & 0xffff);                                    \
    }                                                                          \
    printf("\n");                                                              \
  });
/// Debug tool: print a block in hexa format after izz.
#define print_izz(block)                                                       \
  ({                                                                           \
    printf("\n[   izz] ");                                                     \
    for (uint8_t i = 0; i < BLOCK_LEN; i++) {                                  \
      for (uint8_t j = 0; j < BLOCK_LEN; j++) {                                \
        printf("%x ", block[i][j] & 0xffff);                                   \
      }                                                                        \
    }                                                                          \
    printf("\n");                                                              \
  });
/// Debug tool: print a block in hexa format after idct.
#define print_idct(block)                                                      \
  ({                                                                           \
    printf("\n[  idct] ");                                                     \
    for (uint8_t i = 0; i < BLOCK_LEN; i++) {                                  \
      for (uint8_t j = 0; j < BLOCK_LEN; j++) {                                \
        printf("%x ", block[i][j] & 0xff);                                     \
      }                                                                        \
    }                                                                          \
    printf("\n");                                                              \
  });
/// Debug tool: print a YCbCr mcu in hexa format.
#define print_ycbcr_mcu(mcu, component)                                        \
  ({                                                                           \
    printf("\n* component mcu\n[   mcu] ");                                    \
    for (uint8_t i = 0; i < mcu->height; i++) {                                \
      for (uint8_t j = 0; j < mcu->width; j++) {                               \
        if (component == Y) {                                                  \
          printf("%x ", mcu->blockY[i][j] & 0xff);                             \
        } else if (component == Cb) {                                          \
          printf("%x ", mcu->blockCb[i][j] & 0xff);                            \
        } else {                                                               \
          printf("%x ", mcu->blockCr[i][j] & 0xff);                            \
        }                                                                      \
      }                                                                        \
    }                                                                          \
    printf("\n");                                                              \
  });
#endif
