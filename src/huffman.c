#include "../include/huffman.h"
#include <stdlib.h>

#define MAX_HUFFMAN_DEPTH 16

/// A node in a the Huffman table (tree).
struct huff_table {
  struct huff_table **down;
  bool has_value;
  uint8_t val;
};

/// Builds and returns a empty tree.
struct huff_table *huffman_gen() {
  struct huff_table *table =
      (struct huff_table *)malloc(sizeof(struct huff_table));
  if (table == NULL) {
    abort();
  }

  table->down = (struct huff_table **)malloc(2 * sizeof(struct huff_table *));
  if (table->down == NULL) {
    abort();
  }

  table->down[0] = NULL;
  table->down[1] = NULL;

  table->has_value = false;
  table->val = 0;

  return table;
}

/// Fills the tree.
struct huff_table *get_next_empty_table(struct huff_table **tables,
                                        uint16_t nb_tables) {
  for (uint16_t i = 0; i < nb_tables; i++) {
    for (uint8_t j = 0; j < 2; j++) { // right & left child tables
      if (tables[i]->down[j] == NULL) {
        tables[i]->down[j] = huffman_gen();
        return tables[i]->down[j];
      }
    }
  }
  return NULL;
}

/// Reads the stream and builds a Huffman table. At the beginning, the stream
/// must be 3 bytes after the flag DHT.
struct huff_table *huffman_load_table(struct bitstream *stream,
                                      uint16_t *nb_byte_read) {
  struct huff_table *table = huffman_gen();
  struct huff_table *new_table = NULL;
  struct huff_table *missing_tables = NULL;

  struct huff_table **curr_tables =
      (struct huff_table **)malloc(sizeof(struct huff_table *));
  if (curr_tables == NULL) {
    abort();
  }
  struct huff_table **next_tables = NULL;

  uint8_t values[MAX_HUFFMAN_DEPTH];

  uint32_t bytes;
  uint16_t nb_tables = 1, nb_empty_tables;
  uint8_t computed_max_depth = 0;

  *nb_byte_read = 0;

  for (uint8_t i = 0; i < MAX_HUFFMAN_DEPTH; i++) {
    bitstream_read(stream, 8, &bytes, false);
    values[i] = bytes;
    (*nb_byte_read)++;
  }

  *curr_tables = table;

  for (uint8_t i = 0; i < MAX_HUFFMAN_DEPTH; i++) {
    if (values[i] != 0) {
      computed_max_depth = i;
    }
  }
  computed_max_depth++;

  for (uint8_t i = 0; i < computed_max_depth; i++) {
    for (uint8_t j = 0; j < values[i]; j++) { // copy first values
      bitstream_read(stream, 8, &bytes, false);
      (*nb_byte_read)++;

      new_table = get_next_empty_table(curr_tables, nb_tables);

      if (new_table != NULL) {
        new_table->has_value = true;
        new_table->val = bytes;
      }
    }

    if (i != computed_max_depth - 1) { // shallower nodes
      nb_empty_tables = 2 * nb_tables - values[i];

      next_tables = (struct huff_table **)malloc(nb_empty_tables *
                                                 sizeof(struct huff_table *));
      if (next_tables == NULL) {
        abort();
      }

      missing_tables = get_next_empty_table(curr_tables, nb_tables);

      for (uint8_t k = 0; missing_tables != NULL; k++) {
        next_tables[k] = missing_tables;
        missing_tables = get_next_empty_table(curr_tables, nb_tables);
      }
      free(curr_tables);
      curr_tables = next_tables;
      nb_tables = nb_empty_tables;
    }
  }
  if (computed_max_depth == 1) {
    free(curr_tables);
  }
  if (next_tables != NULL) {
    free(next_tables);
  }
  return table;
}

/// Returns the next value reached by using the Huffman table according to the
/// bits read in the stream.
int8_t huffman_next_value(struct huff_table *table, struct bitstream *stream) {
  uint32_t byte = 0;
  while (!table->has_value) {
    bitstream_read(stream, 1, &byte, true);
    table = table->down[byte];
  }
  return table->val;
}

/// Frees the memory used by the Huffman tables.
void huffman_free_table(struct huff_table *table) {
  if (table != NULL) {
    huffman_free_table(table->down[0]);
    huffman_free_table(table->down[1]);
    free(table->down);
    free(table);
  }
}
