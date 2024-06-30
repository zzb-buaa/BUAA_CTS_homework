#ifndef ALGORITHM_H_
#define ALGORITHM_H_
#define true 1
#define false 0
#define FS 100
#define BUFFER_SIZE  (FS * 5)
#define HR_FIFO_SIZE 7
#define MA4_SIZE  4 // DO NOT CHANGE
#define HAMMING_SIZE  5 // DO NOT CHANGE
#define min(x,y) ((x) < (y) ? (x) : (y))
#include <stdint.h>

extern const uint16_t auw_hamm[31]; // Hamm=  long16(512* hamming(5)');
extern const uint8_t uch_spo2_table[184];
static int32_t an_dx[BUFFER_SIZE - MA4_SIZE]; // delta
static int32_t an_x[BUFFER_SIZE]; // ir

void maxim_heart_rate_only(uint32_t *pun_ir_buffer, int32_t n_ir_buffer_length, int32_t *pn_heart_rate, int8_t *pch_hr_valid);
void maxim_find_peaks(int32_t *pn_locs, int32_t *pn_npks, int32_t *pn_x, int32_t n_size, int32_t n_min_height, int32_t n_min_distance, int32_t n_max_num);
void maxim_peaks_above_min_height(int32_t *pn_locs, int32_t *pn_npks, int32_t *pn_x, int32_t n_size, int32_t n_min_height);
void maxim_remove_close_peaks(int32_t *pn_locs, int32_t *pn_npks, int32_t *pn_x, int32_t n_min_distance);
void maxim_sort_ascend(int32_t *pn_x, int32_t n_size);
void maxim_sort_indices_descend(int32_t *pn_x, int32_t *pn_indx, int32_t n_size);

#endif /* ALGORITHM_H_ */

