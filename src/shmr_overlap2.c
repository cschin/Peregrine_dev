#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <wordexp.h>
#include "kalloc.h"
#include "khash.h"
#include "kvec.h"
#include "shimmer.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define handle_error(msg) \
  do {                    \
    perror(msg);          \
    exit(EXIT_FAILURE);   \
  } while (0)

#define MMER_COUNT_LOWER_BOUND 2
#define MMER_COUNT_UPPER_BOUND 4096
#ifndef ORIGINAL
#define ORIGINAL 0
#endif
#ifndef REVERSED
#define REVERSED 1
#endif
#define READ_END_FUZZINESS 48
#define LOCAL_OVERLAP_UPPERBOUND 0
#define BESTN 5
#define OVERLAP 0
#define CONTAINS 1
#define CONTAINED 2
#define ALNBANDSIZE 100

ovlp_match_t *ovlp_match2(uint8_t *query_seq, seq_coor_t q_len, uint8_t q_strand,
                         uint8_t *target_seq, seq_coor_t t_len,
                         uint8_t t_strand, seq_coor_t band_tolerance,
                         bool get_deltas,
                         bool use_hpc);

typedef struct {
    uint32_t rid0, rid1, strand1, len0;
    int32_t d_left, d_right;
} ovlp_candidate_t;


typedef struct {
  size_t n, m;
  ovlp_candidate_t *a;
} ovlp_candidate_v;

void build_map2(mm128_v *mmers, khash_t(MMER0) * mmer0_map,
               khash_t(RLEN) * rlmap, khash_t(MMC) * mcmap, uint32_t chunk,
               uint32_t total_chunk, uint32_t lowerbound, uint32_t upperbound) {
  uint64_t mhash;
  mm128_t mmer0, mmer1;
  mp128_v *mpv;
  uint32_t rid;
  uint32_t pos, rpos;
  uint32_t span;
  uint32_t mcount = 0;
  int32_t absent;
  mp128_t rp;
  khiter_t k;
  khash_t(MMER1) * mmer1_map;
  size_t s = 0;

  for (;;) {
    if (s >= mmers->n) break;
    mmer0 = mmers->a[s];
    mhash = mmer0.x >> 8;
    k = kh_get(MMC, mcmap, mhash);
    assert(k != kh_end(mcmap));
    mcount = kh_val(mcmap, k);
    if (mcount >= lowerbound && mcount < upperbound) break;
    s++;
  }

  for (size_t i = s + 1; i < mmers->n; i++) {
    mmer1 = mmers->a[i];
    mhash = mmer1.x >> 8;
    k = kh_get(MMC, mcmap, mhash);
    assert(k != kh_end(mcmap));
    mcount = kh_val(mcmap, k);
    if (mcount < lowerbound || mcount > upperbound) continue;

    if ((mmer0.y >> 32) != (mmer1.y >> 32)) {  // the pairs are not in the same read
       mmer0 = mmer1;
       continue;
    }

    if ((mmer0.y >> 32) % total_chunk != chunk % total_chunk) {
      mmer0 = mmer1;
      continue;
    }

    // don't use two minimers that are too close to each other
    if (((mmer1.y >> 1) & 0xFFFFFFF) - ((mmer0.y >> 1) & 0xFFFFFFF) < 100) {
      mmer0 = mmer1;
      continue;
    }

    //DEBUG
    //if ((mmer1.y >> 32) == 2130) {
    //   printf("DEBUG mmcount: %d %d %d\n", mmer1.y >> 32, ((mmer1.y >> 1) & 0xFFFFFFF), mcount);
    //}
    //DEBUG END

    // build the index for reads in the chunk
    k = kh_put(MMER0, mmer0_map, mmer0.x, &absent);
    if (absent) kh_value(mmer0_map, k) = kh_init(MMER1);

    mmer1_map = kh_value(mmer0_map, k);
    k = kh_put(MMER1, mmer1_map, mmer1.x, &absent);
    if (absent) {
      mpv = kmalloc(NULL, sizeof(mp128_v));
      mpv->n = 0;
      mpv->m = 0;
      mpv->a = NULL;
      kh_value(mmer1_map, k) = mpv;
    } else {
      mpv = kh_value(mmer1_map, k);
    }
    rp.y0 = mmer0.y;
    rp.y1 = mmer1.y;

    rp.direction = ORIGINAL;
    kv_push(mp128_t, NULL, *mpv, rp);
    mmer0 = mmer1;
  }
}

mp128_v * match_shimmer_pair(khash_t(MMER0) * mmer0_map, 
                             mm128_t mmer0,
                             mm128_t mmer1) {
  khiter_t k;
  mp128_v *mpv;
  khash_t(MMER1) * mmer1_map;
  k = kh_get(MMER0, mmer0_map, mmer0.x);
  if (k == kh_end(mmer0_map)) {
      return NULL;
  }
  mmer1_map = kh_value(mmer0_map, k);
  k = kh_get(MMER1, mmer1_map, mmer1.x);
  if (k == kh_end(mmer1_map)) {
      return NULL;
  }
  mpv = kh_value(mmer1_map, k);
  return mpv;
}

KHASH_MAP_INIT_INT64(RPAIR, uint32_t);

uint32_t get_rid_pair_count(khash_t(RPAIR) *rid_pairs, 
                    uint32_t rid0,
                    uint32_t rid1) {
  uint64_t ridp;
  khiter_t k;
  ridp = (((uint64_t)rid0) << 32) | ((uint64_t)rid1);
  k = kh_get(RPAIR, rid_pairs, ridp);
  if (k == kh_end(rid_pairs)) {
      return 0;
  } else {
      return kh_val(rid_pairs, k);
  }
}

void increase_rid_pair_count(khash_t(RPAIR) *rid_pairs,
                  uint32_t rid0,
                  uint32_t strand0,
                  uint32_t rid1,
                  uint32_t strand1) {
  uint64_t ridp;
  int32_t absent;
  khiter_t k;
  uint32_t r0 = (rid0 << 1) | strand0;
  uint32_t r1 = (rid1 << 1) | strand1;

  ridp = (((uint64_t)r0) << 32) | ((uint64_t)r1);
  k = kh_get(RPAIR, rid_pairs, ridp);
  if (k == kh_end(rid_pairs)) {
     k = kh_put(RPAIR, rid_pairs, ridp, &absent);
     kh_val(rid_pairs, k) = 1;
  } else {
     kh_val(rid_pairs, k) += 1;
  }
}

// get precise d_left from an estimation
#define ALNW 200
#define ALNH 24
int32_t update_d_left(uint8_t * seq0, uint8_t * seq1, 
                          int32_t d_left, 
                          uint32_t rlen0, uint32_t rlen1,
                          uint32_t strand1) {
    int32_t score[ALNW+1][ALNH+1];
    memset(score, 0, (ALNW+1)*(ALNH+1)*sizeof(int32_t));
    uint32_t xs, xe;
    if (d_left >= 0) {
      xs = d_left - (ALNW>>1) > 0 ? d_left - (ALNW>>1) : 0;
      xe = d_left + (ALNW>>1) < rlen0 ? d_left + (ALNW>>1) : rlen0;
    } else {
      xs = abs(d_left) - (ALNW>>1) > 0 ? abs(d_left) - (ALNW>>1) : 0;
      xe = abs(d_left) + (ALNW>>1) < rlen1 ? abs(d_left) + (ALNW>>1) : rlen1;
    }
    uint32_t max_x = 0;
    uint32_t max_y = 0;
    int32_t max_s = 0;
    for (int32_t y = ALNH - 1 ; y >= 0; y--) {
      for (int32_t x = xe - xs - 1; x >= 0; x--) {
        uint8_t b0, b1;  
        if (d_left > 0) {
          b0 = seq0[x+xs] & 0x0F;
          b1 = strand1 == 0 ? seq1[y] : seq1[y] >> 4;
          b1 &= 0x0F; 
        } else {
          b0 = seq0[y] & 0x0F;
          b1 = strand1 == 0 ? seq1[x+xs] : seq1[x+xs] >> 4;
          b1 &= 0x0F; 
        }
            
        int32_t s0, sd, sy, sx;
        s0 = b0 == b1 ? 2 : -4;
        sy = score[x][y+1] - 4;
        sd = s0 + score[x+1][y+1];
        sx = score[x+1][y] - 4; 
        if (sd > sx && sd > sy) {
           score[x][y] = sd;
        } else if (sx > sy) {
           score[x][y] = sx;
        } else {
           score[x][y] = sy;
        }
        if (score[x][y] > max_s) {
          max_s = score[x][y];
          max_x = x;
          max_y = y; 
        } 
      }
   }
  if (max_y < 2 && max_s >= (ALNH-2) * 2) {
      if (d_left > 0) {
          return xs + max_x;
      } else {
          return -xs-max_x;
      }
  } else {
      return 0;
  }
}

ovlp_match_t * match_seqs(uint8_t * seq0, uint8_t * seq1, 
                          int32_t d_left, 
                          uint32_t rlen0, uint32_t rlen1,
                          uint32_t strand1) {
      //printf("DEBUG: d_left:%d rlen1:%d rlen1:%d\n", d_left, rlen0, rlen1);
      assert(d_left < (int32_t) rlen0);
      assert(d_left + (int32_t) rlen1 > 0);
      ovlp_match_t *match;
      uint32_t align_bandwidth = 100;
      if (d_left > 0 ) {
        match = ovlp_match2(seq0 + d_left, rlen0 - d_left, ORIGINAL, 
                           seq1, rlen1, strand1, 
                           align_bandwidth, true, true);
      } else {
        d_left = abs(d_left);
        match = ovlp_match2(seq0, rlen0, ORIGINAL, 
                           seq1 + d_left, rlen1 - d_left, strand1, 
                           align_bandwidth, true,  true);
      }
      return match;
}

KHASH_MAP_INIT_INT(OVLP_CANDIDATES, ovlp_candidate_v *);

void push_ovlp_candidate(khash_t(OVLP_CANDIDATES) * ovlp_candidates,
                         ovlp_candidate_t *candidate) {
  uint32_t rid0 = candidate->rid0;
  int32_t absent;
  khiter_t k;
  //DEBUG
  //if (rid0 == 985) {
  //    printf("DEBUG X: %d %d %d \n", rid0, candidate->rid1, candidate->d_left);
  //}
  //DEBUG END
  k = kh_get(OVLP_CANDIDATES, ovlp_candidates, rid0);
  if (k == kh_end(ovlp_candidates)) {
     k = kh_put(OVLP_CANDIDATES, ovlp_candidates, rid0, &absent);
     ovlp_candidate_v *v;
     v = (ovlp_candidate_v *) calloc(1, sizeof(ovlp_candidate_v));
     kv_init(*v);
     kh_val(ovlp_candidates, k) = v;
     kv_push(ovlp_candidate_t, NULL, *v, *candidate);
  } else {
     ovlp_candidate_v *v;
     v = kh_val(ovlp_candidates, k);
     kv_push(ovlp_candidate_t, NULL, *v, *candidate);
  }
}

int ovlp_comp_left(const void *a, const void *b) {
  ovlp_candidate_t *a0 = (ovlp_candidate_t *)a;
  ovlp_candidate_t *b0 = (ovlp_candidate_t *)b;
  return (int32_t) a0->d_left >= (int32_t) b0->d_left;
}

int ovlp_comp_right(const void *a, const void *b) {
  ovlp_candidate_t *a0 = (ovlp_candidate_t *)a;
  ovlp_candidate_t *b0 = (ovlp_candidate_t *)b;
  return (int32_t) a0->d_right <= (int32_t) b0->d_right;
}

int ovlp_comp_rid(const void *a, const void *b) {
  ovlp_candidate_t *a0 = (ovlp_candidate_t *)a;
  ovlp_candidate_t *b0 = (ovlp_candidate_t *)b;
  return (a0->rid1) == (b0->rid1) ? (int32_t) a0->d_left >= (int32_t) b0->d_left : (a0->rid1) > (b0->rid1);
}

bool check_match(ovlp_match_t *match, uint32_t slen0, uint32_t slen1) {
  seq_coor_t q_bgn, q_end, t_bgn, t_end;
  q_bgn = match->q_bgn;
  q_end = match->q_end;
  t_bgn = match->t_bgn;
  t_end = match->t_end;
  double err_est;
  if (match->m_size == 0) return false;
  err_est = 100.0 - 100.0 * (double)(match->dist) / (double)(match->m_size+1);
  // printf("DEBUG: check_match: %d %d %d %d %d %d %f\n", q_bgn, q_end, slen0, t_bgn, t_end, slen1, err_est);
  if (err_est < 99.5) return false;;
  if (q_end < 500 || t_end < 500)  return false;

  if (q_bgn > READ_END_FUZZINESS || t_bgn > READ_END_FUZZINESS) return false; 
  if ((abs((int64_t) (slen0) - (int64_t) q_end) > READ_END_FUZZINESS) &&  
      (abs((int64_t) (slen1) - (int64_t) t_end) > READ_END_FUZZINESS)) return false;
  return true;
}

KHASH_MAP_INIT_INT64(MRID, uint32_t);

void dump_candidates(khash_t(OVLP_CANDIDATES) * ovlp_candidates,
                     khash_t(RLEN) * rlmap,
                     uint32_t ovlp_upper,
                     uint32_t bestn,
                     uint8_t *seq_p,
                     FILE * ovlp_file) {
  ovlp_candidate_v * v;
  uint32_t rlen0,  rlen1;
  uint32_t rid0;
  uint8_t *seq0 = NULL;
  uint8_t *seq1 = NULL;
  khiter_t k;
  ovlp_match_t *match;

  for (khiter_t __i = kh_begin(ovlp_candidates); 
         __i != kh_end(ovlp_candidates); 
         ++__i) {

    khash_t(MRID) *mrid = kh_init(MRID);
    int32_t absent;

    if (!kh_exist(ovlp_candidates, __i)) continue;
    v = kh_val(ovlp_candidates, __i);

    if (ovlp_upper != 0) {
       if (v->n > ovlp_upper) continue;
    }

    rid0 = v->a[0].rid0;
    k = kh_get(RLEN, rlmap, rid0);
    assert(k != kh_end(rlmap));
    rlen0 = kh_val(rlmap, k).len;
    seq0 = get_read_seq_mmap_ptr(seq_p, rid0, rlmap);
    qsort(v->a, v->n, sizeof(ovlp_candidate_t), ovlp_comp_left);
    uint32_t ovlp_count=0;
    int32_t last_d_left = -1;
    uint32_t current_rid1 = v->a[0].rid1;
    for (uint32_t i = 0; i < v->n; i++) {
      ovlp_candidate_t c;
      c = v->a[i];
      if (c.d_left <= 0) continue;

      if (c.rid1 != current_rid1) {
          last_d_left = -1;
          current_rid1 = c.rid1;
      }


      k = kh_get(RLEN, rlmap, c.rid1);
      assert(k != kh_end(rlmap));
      rlen1 = kh_val(rlmap, k).len;

      seq1 = get_read_seq_mmap_ptr(seq_p, c.rid1, rlmap);

      uint32_t tmp;
      tmp = c.d_left;
      c.d_left = update_d_left(seq0, seq1, 
                         c.d_left,
                         rlen0, rlen1, 
                         c.strand1);

      if (c.d_left == 0) continue;

      k = kh_get(MRID, mrid, (((uint64_t) c.rid1) << 32 | ((uint64_t) c.d_left)));
      if (k != kh_end(mrid)) continue;

      c.d_right += c.d_left - tmp;

      //if (abs(c.d_left - last_d_left) < 24) {
      //    last_d_left = c.d_left;
      //    continue;
      //}

      //DEBUG
      // printf("DEBUG:1 rid0, rid1= %d %d %d\n", c.rid0, c.rid1, c.d_left);
      //DEBUG END
      
      last_d_left = c.d_left;
      

      match = match_seqs(seq0, seq1, 
                         c.d_left,
                         rlen0, rlen1, 
                         c.strand1);
      if (check_match(match, rlen0 - c.d_left, rlen1) == false) {
          k = kh_put(MRID, mrid, (((uint64_t) c.rid1) << 32 | ((uint64_t) c.d_left)), &absent);
          kh_val(mrid, k) = 1;
          free_ovlp_match(match);
          continue;
      } 
      seq_coor_t q_bgn, q_end, t_bgn, t_end;
      q_bgn = match->q_bgn;
      q_end = match->q_end;
      t_bgn = match->t_bgn;
      t_end = match->t_end;

      if (c.d_right > 2500) ovlp_count++;
      double err_est;
      err_est = 100.0 - 100.0 * (double)(match->dist) / (double)(match->m_size);
      
      //DEBUG
      /*
      printf("DEBUG:1 matched rid0, rid1= %d %d %d\n", c.rid0, c.rid1, c.d_left);
      printf( "%d %d %d %d %d %d %d %d %d %d %0.2f\n",
             c.rid0,
             c.rid1,
             c.strand1,
             c.len0,
             c.d_left,
             c.d_right,
             c.d_left + q_bgn,
             c.d_left + q_end,
             t_bgn,
             t_end,
             err_est);
      printf("\n");
      */
      //DEBUG END
      
      fprintf(ovlp_file,
             "O %d %d %d %d %d %d %d %d %d %d %d %0.2f\n",
             c.rid0,
             c.rid1,
             c.strand1,
             c.len0,
             c.d_left,
             c.d_right,
             c.d_left + q_bgn,
             c.d_left + q_end,
             t_bgn,
             t_end,
             match->dist,
             err_est);
      fprintf(ovlp_file, "D %d %d %d", c.rid0, c.rid1,  match->dist);
      for (uint32_t idx = 0; idx < match->dist; idx++) {
        uint32_t val = match->reduce_deltas[idx];
        fprintf(ovlp_file, " %d %d", c.d_left + (val>>1), val&0x1);
      }
      fprintf(ovlp_file, "\n");
      free_ovlp_match(match);

      k = kh_put(MRID, mrid, (((uint64_t) c.rid1) << 32 | ((uint64_t) c.d_left)), &absent);
      kh_val(mrid, k) = 1;

      if (ovlp_count > bestn) break;
    }
    qsort(v->a, v->n, sizeof(ovlp_candidate_t), ovlp_comp_right);
    ovlp_count=0;
    int32_t last_d_right = 1;
    current_rid1 = v->a[0].rid1;
    for (uint32_t i = 0; i < v->n; i++) {
      ovlp_candidate_t c;
      c = v->a[i];
      if (c.d_right >= 0) continue;
      
      if (c.rid1 != current_rid1) {
          last_d_right = 1;
          current_rid1 = c.rid1;
      }

      k = kh_get(RLEN, rlmap, c.rid1);
      assert(k != kh_end(rlmap));
      rlen1 = kh_val(rlmap, k).len;

      seq1 = get_read_seq_mmap_ptr(seq_p, c.rid1, rlmap);
      uint32_t tmp;
      tmp = c.d_left;
      c.d_left = update_d_left(seq0, seq1, 
                         c.d_left,
                         rlen0, rlen1, 
                         c.strand1);

      if (c.d_left == 0) continue;

      k = kh_get(MRID, mrid, (((uint64_t) c.rid1) << 32 | ((uint64_t) c.d_left)));
      if (k != kh_end(mrid)) continue;

      c.d_right += c.d_left - tmp;

      //if (abs(c.d_right - last_d_right) < 24) {
      //    last_d_right = c.d_right;
      //    continue;
      //}
      //DEBUG
      //printf("DEBUG:2 rid0, rid1= %d %d %d\n", c.rid0, c.rid1, c.d_right);
      //DEBUG END

      last_d_right = c.d_right;
      
      match = match_seqs(seq0, seq1, 
                         c.d_left,
                         rlen0, rlen1, 
                         c.strand1);
      //DEBUG
      //if (c.rid0 == 1339 && c.rid1 == 1217) {
      //  printf("DEBUG:3 rid0, rid1= %d %d\n", c.rid0, c.rid1);
      //  printf("d_left, rlen0, rlen1, c.strnd1= %d %d %d %d\n", c.d_left, rlen0, rlen1, c.strand1);
      //  printf("%d %d %d %d\n", q_bgn, q_end, t_bgn, t_end);
      //}
      //DEBUG END

      if (check_match(match, rlen0, rlen1 - abs(c.d_left)) == false) {
          k = kh_put(MRID, mrid, (((uint64_t) c.rid1) << 32 | ((uint64_t) c.d_left)), &absent);
          kh_val(mrid, k) = 1;
          free_ovlp_match(match);
          continue;
      }       
      seq_coor_t q_bgn, q_end, t_bgn, t_end;
      q_bgn = match->q_bgn;
      q_end = match->q_end;
      t_bgn = match->t_bgn;
      t_end = match->t_end;
      if (c.d_left < -2500) ovlp_count++;
      double err_est;
      err_est = 100.0 - 100.0 * (double)(match->dist) / (double)(match->m_size);
      fprintf(ovlp_file,
              "O %d %d %d %d %d %d %d %d %d %d %d %0.2f\n",
              c.rid0,
              c.rid1,
              c.strand1,
              c.len0,
              c.d_left,
              c.d_right,
              q_bgn,
              q_end,
              abs(c.d_left) + t_bgn,
              abs(c.d_left) + t_end,
              match->dist,
              err_est);
      fprintf(ovlp_file, "D %d %d %d", c.rid0, c.rid1,  match->dist);
      for (uint32_t idx = 0; idx < match->dist; idx++) {
        uint32_t val = match->reduce_deltas[idx];
        fprintf(ovlp_file, " %d %d", val>>1, val&0x1);
      }
      fprintf(ovlp_file, "\n");
      free_ovlp_match(match);

      k = kh_put(MRID, mrid, (((uint64_t) c.rid1) << 32 | ((uint64_t) c.d_left)), &absent);
      kh_val(mrid, k) = 1;

      if (ovlp_count > bestn) break;
    }
    kh_destroy(MRID, mrid);
  }
}

void filter_candidates(khash_t(OVLP_CANDIDATES) * in_ovlps,
        khash_t(OVLP_CANDIDATES) * out_ovlps) {
  for (khiter_t __i = kh_begin(in_ovlps); 
         __i != kh_end(in_ovlps); 
         ++__i) {

    khash_t(MRID) *mrid = kh_init(MRID);
    int32_t absent;
    ovlp_candidate_v * v;

    if (!kh_exist(in_ovlps, __i)) continue;
    v = kh_val(in_ovlps, __i);
    qsort(v->a, v->n, sizeof(ovlp_candidate_t), ovlp_comp_rid);
    uint32_t ovlp_count=0;
    int32_t last_d_left = -1;
    uint32_t current_rid1 = v->a[0].rid1;
    uint32_t cluster_count = 0;
    ovlp_candidate_t c;
    ovlp_candidate_t last_c;
    for (uint32_t i = 0; i < v->n; i++) {
      c = v->a[i];
      if (last_d_left == -1 || current_rid1 != c.rid1) {
          last_d_left = c.d_left;
          if (cluster_count > 3) {
            push_ovlp_candidate(out_ovlps, &last_c);
          }
          cluster_count = 0;
          current_rid1 = c.rid1;
          last_c = c;
      }
      if (c.d_left < last_d_left + 72) {
        cluster_count++;
      } else {
        if (cluster_count > 3) {
          push_ovlp_candidate(out_ovlps, &last_c);
        }
        cluster_count = 1;
        last_d_left = c.d_left;
      }
      //DEBUG
      /*
      printf( "%d %d %d %d %d %d %d\n",
              c.rid0,
              c.rid1,
              c.strand1,
              c.len0,
              c.d_left,
              c.d_right,
              cluster_count);
      */
      last_c = c;
    }
    if (cluster_count > 3) {
      push_ovlp_candidate(out_ovlps, &last_c);
    }
  }
}


void build_ovlp_candidates(mm128_v *mmers, 
               khash_t(MMER0) * mmer0_map,
               khash_t(RLEN) * rlmap, 
               khash_t(MMC) * mcmap, 
               uint32_t chunk,
               uint32_t total_chunk, 
               uint32_t lowerbound, 
               uint32_t upperbound,
               khash_t(OVLP_CANDIDATES) * ovlp_candidates) {
  uint64_t mhash;
  mm128_t mmer0, mmer1;
  mp128_v *mpv;
  uint32_t rid;
  uint32_t pos, rpos;
  uint32_t span;
  uint32_t mcount = 0;
  int32_t absent;
  mp128_t rp;
  khiter_t k;
  khash_t(MMER1) * mmer1_map;
  size_t s = 0;
  // uint8_t *seq0 = NULL;
  // uint8_t *seq1 = NULL;
  ovlp_match_t *match;
  ovlp_candidate_t candidate;
  khash_t(OVLP_CANDIDATES) * ovlp_tmp = kh_init(OVLP_CANDIDATES);
  khash_t(RPAIR) *rid_pairs = kh_init(RPAIR);

  for (;;) {
    if (s >= mmers->n) break;
    mmer0 = mmers->a[s];
    mhash = mmer0.x >> 8;
    k = kh_get(MMC, mcmap, mhash);
    assert(k != kh_end(mcmap));
    mcount = kh_val(mcmap, k);
    if (mcount >= lowerbound && mcount < upperbound) break;
    s++;
  }

  for (size_t i = s + 1; i < mmers->n; i++) {
    uint32_t rid0, pos0, rlen0, strand0;
    uint32_t rid1, pos1, rlen1, strand1;

    mmer1 = mmers->a[i];
    mhash = mmer1.x >> 8;
    k = kh_get(MMC, mcmap, mhash);
    assert(k != kh_end(mcmap));
    mcount = kh_val(mcmap, k);
    if (mcount < lowerbound || mcount > upperbound) continue;

    if ((mmer0.y >> 32) != (mmer1.y >> 32)) {  // the pairs are in the same read
      filter_candidates(ovlp_tmp, ovlp_candidates);
      for (khiter_t __i = kh_begin(ovlp_tmp); 
           __i != kh_end(ovlp_tmp); 
           ++__i) {
        ovlp_candidate_v * v;
        if (!kh_exist(ovlp_tmp, __i)) continue;
        v = kh_val(ovlp_tmp, __i);
        kv_destroy(*v);
      }
      kh_destroy(OVLP_CANDIDATES, ovlp_tmp);
      ovlp_tmp = kh_init(OVLP_CANDIDATES);

      kh_destroy(RPAIR, rid_pairs);
      rid_pairs = kh_init(RPAIR);
      mmer0 = mmer1;
      continue;
    }

      // don't use two minimers that are too close to each other
    if (((mmer1.y >> 1) & 0xFFFFFFF) - ((mmer0.y >> 1) & 0xFFFFFFF) < 100) {
      mmer0 = mmer1;
      continue;
    }

    mpv = match_shimmer_pair(mmer0_map, mmer0, mmer1);
    if (mpv != NULL) {
      rid1 = (uint32_t)(mmer0.y >> 32);
      pos1 = (uint32_t)((mmer0.y & 0xFFFFFFFF) >> 1) + 1;
      k = kh_get(RLEN, rlmap, rid1);
      assert(k != kh_end(rlmap));
      rlen1 = kh_val(rlmap, k).len;
      strand1 = ORIGINAL;
      for (size_t __k0 = 0; __k0 <  (mpv->n) ; __k0++) { 
        uint64_t y0;
        y0 = mpv->a[__k0].y0;
        rid0 = (uint32_t)(y0 >> 32);
        if (rid0 == rid1) {
            continue;
        }
        if (get_rid_pair_count(rid_pairs, rid0, rid1) > 12) continue;

        k = kh_get(RLEN, rlmap, rid0);
        assert(k != kh_end(rlmap));
        rlen0 = kh_val(rlmap, k).len;

        pos0 = (uint32_t)((y0 & 0xFFFFFFFF) >> 1) + 1;

        candidate.rid0 = rid0;
        candidate.rid1 = rid1;
        candidate.strand1 = strand1;
        candidate.len0 = rlen0;
        candidate.d_left = (int32_t) pos0 - (int32_t) pos1;
        candidate.d_right = (int32_t) pos0 - (int32_t) pos1 + (int32_t) rlen1 - (int32_t) rlen0;
        push_ovlp_candidate(ovlp_tmp, &candidate);
        increase_rid_pair_count(rid_pairs, rid0, 0, rid1, strand1);
      }
    }

    // revers
    mpv = match_shimmer_pair(mmer0_map, mmer1, mmer0);
    if (mpv != NULL) {
      rid1 = (uint32_t)(mmer1.y >> 32);
      k = kh_get(RLEN, rlmap, rid1);
      assert(k != kh_end(rlmap));
      rlen1 = kh_val(rlmap, k).len;
      pos1 = (uint32_t)((mmer1.y & 0xFFFFFFFF) >> 1) + 1;
      span = mmer1.x & 0xFF;
      pos1 = rlen1 - pos1 + span - 1;
      strand1 = REVERSED;
      for (size_t __k0 = 0; __k0 <  (mpv->n) ; __k0++) { 
        uint64_t y0;
        y0 = mpv->a[__k0].y0;
        rid0 = (uint32_t)(y0 >> 32);
        if (rid0 == rid1) {
          continue;
        }
        if (get_rid_pair_count(rid_pairs, rid0, rid1) > 12) continue;

        k = kh_get(RLEN, rlmap, rid0);
        assert(k != kh_end(rlmap));
        rlen0 = kh_val(rlmap, k).len;

        pos0 = (uint32_t)((y0 & 0xFFFFFFFF) >> 1) + 1;
        
        candidate.rid0 = rid0;
        candidate.rid1 = rid1;
        candidate.strand1 = strand1;
        candidate.len0 = rlen0;
        candidate.d_left = (int32_t) pos0 - (int32_t) pos1;
        candidate.d_right = (int32_t) pos0 - (int32_t) pos1 + (int32_t) rlen1 - (int32_t) rlen0;
        push_ovlp_candidate(ovlp_tmp, &candidate);
        increase_rid_pair_count(rid_pairs, rid0, 0, rid1, strand1);
      }
    }
    mmer0 = mmer1;
  }

  filter_candidates(ovlp_tmp, ovlp_candidates);

  for (khiter_t __i = kh_begin(ovlp_tmp); 
         __i != kh_end(ovlp_tmp); 
         ++__i) {
    ovlp_candidate_v * v;
    if (!kh_exist(ovlp_tmp, __i)) continue;
    v = kh_val(ovlp_tmp, __i);
    kv_destroy(*v);
  }
  kh_destroy(OVLP_CANDIDATES, ovlp_tmp);
  kh_destroy(RPAIR, rid_pairs);
}


int mp128_comp(const void *a, const void *b) {
  mp128_t *a0 = (mp128_t *)a;
  mp128_t *b0 = (mp128_t *)b;
  return ((a0->y0 & 0xFFFFFFFF) >> 1) < ((b0->y0 & 0xFFFFFFFF) >> 1);
}


int main(int argc, char *argv[]) {
  char *seqdb_prefix = NULL;
  char *shimmer_prefix = NULL;
  char *ovlp_file_path = NULL;

  FILE *ovlp_file;

  char mmc_file_path[8192];
  char mmer_file_path[8192];
  char seq_idx_file_path[8192];
  char seqdb_file_path[8291];
  int c;
  uint32_t bestn = BESTN;
  uint32_t mc_lower = MMER_COUNT_LOWER_BOUND;
  uint32_t mc_upper = MMER_COUNT_UPPER_BOUND;
  uint32_t ovlp_upper = LOCAL_OVERLAP_UPPERBOUND;
  uint32_t align_bandwidth = ALNBANDSIZE;

  uint32_t total_chunk = 1, mychunk = 1;

  wordexp_t p;
  char **mmc_fns;
  char **shimmer_fns;

  mm128_v mmers = {0, 0, 0};
  mm128_v mmers_;
  mm_count_v mmc;

  khash_t(RLEN) * rlmap;
  khash_t(MMC) *mcmap = kh_init(MMC);

  khash_t(MMER0) * mmer0_map;
  khash_t(MMER1) * mmer1_map;

  mp128_v *mpv;

  opterr = 0;

  while ((c = getopt(argc, argv, "p:l:t:c:b:o:m:M:w:n:")) != -1) {
    switch (c) {
      case 'p':
        seqdb_prefix = optarg;
        break;
      case 'l':
        shimmer_prefix = optarg;
        break;
      case 't':
        total_chunk = atoi(optarg);
        break;
      case 'c':
        mychunk = atoi(optarg);
        break;
      case 'b':
        bestn = atoi(optarg);
        break;
      case 'o':
        ovlp_file_path = optarg;
        break;
      case 'm':
        mc_lower = atoi(optarg);
        break;
      case 'M':
        mc_upper = atoi(optarg);
        break;
      case 'w':
        align_bandwidth = atoi(optarg);
        break;
      case 'n':
        ovlp_upper = atoi(optarg);
        break;
      case '?':
        if (optopt == 'p') {
          fprintf(stderr,
                  "Option -%c not specified, using 'seq_dataset' as the "
                  "sequence db prefix\n",
                  optopt);
        }
        if (optopt == 'l') {
          fprintf(stderr,
                  "Option -%c not specified, using 'shimmer-L2' as the L2 "
                  "index prefix\n",
                  optopt);
        }
        if (optopt == 'o') {
          fprintf(
              stderr,
              "Option -%c not specified, using 'ovlp.# as default output' \n",
              optopt);
        }
        return 1;
      default:
        abort();
    }
  }

  assert(total_chunk > 0);
  assert(mychunk > 0 && mychunk <= total_chunk);

  if (seqdb_prefix == NULL) {
    seqdb_prefix = (char *)calloc(8192, 1);
    snprintf(seqdb_prefix, 8191, "seq_dataset");
  }

  if (shimmer_prefix == NULL) {
    shimmer_prefix = (char *)calloc(8192, 1);
    snprintf(shimmer_prefix, 8191, "shimmer-L2");
  }

  if (ovlp_file_path == NULL) {
    ovlp_file_path = (char *)calloc(8192, 1);
    snprintf(ovlp_file_path, 8191, "ovlp.%02d", mychunk);
  }

  int written;
  written = snprintf(seq_idx_file_path, sizeof(seq_idx_file_path), "%s.idx",
                     seqdb_prefix);
  assert(written < sizeof(seq_idx_file_path));
  fprintf(stderr, "using index file: %s\n", seq_idx_file_path);

  rlmap = get_read_length_map(seq_idx_file_path);

  written = snprintf(seqdb_file_path, sizeof(seqdb_file_path), "%s.seqdb",
                     seqdb_prefix);
  assert(written < sizeof(seqdb_file_path));
  fprintf(stderr, "using seqdb file: %s\n", seqdb_file_path);

  written = snprintf(mmer_file_path, sizeof(mmer_file_path),
                     "%s-[0-9]*-of-[0-9]*.dat", shimmer_prefix);
  assert(written < sizeof(mmer_file_path));
  wordexp(mmer_file_path, &p, 0);
  shimmer_fns = p.we_wordv;
  for (int i = 0; i < p.we_wordc; i++) {
    fprintf(stderr, "using shimmer data file: %s\n", shimmer_fns[i]);
    mmers_ = read_mmlist(shimmer_fns[i]);
    append_mmlist(&mmers, &mmers_);
    kv_destroy(mmers_);
  }
  wordfree(&p);

  written = snprintf(mmc_file_path, sizeof(mmc_file_path),
                     "%s-MC-[0-9]*-of-[0-9]*.dat", shimmer_prefix);
  assert(written < sizeof(mmc_file_path));
  wordexp(mmc_file_path, &p, 0);
  mmc_fns = p.we_wordv;

  uint32_t *mc_bin;
  mc_bin = calloc(512, sizeof(uint32_t));
  for (int i = 0; i < p.we_wordc; i++) {
    fprintf(stderr, "using shimmer count file: %s\n", mmc_fns[i]);
    mmc = read_mm_count(mmc_fns[i]);
    for (uint32_t j=0; j < mmc.n; j++) {
      if (mmc.a[j].count > 511) {
          mc_bin[511] += 1;
      } else {
          mc_bin[mmc.a[j].count-1]++;
      }
    }

    aggregate_mm_count(mcmap, &mmc);
    kv_destroy(mmc);
  }
  wordfree(&p);
  
  fprintf(stderr, "shimmer count histogram\n");
  for (uint32_t j=0; j < 512; j++) {
    fprintf(stderr, "%10d %d", j+1, mc_bin[j]);
    fprintf(stderr, "\n");
  }
  free(mc_bin);

  ovlp_file = fopen(ovlp_file_path, "w");

  char buffer[32768];

  setvbuf(ovlp_file, buffer, _IOFBF, sizeof(buffer));

  mmer0_map = kh_init(MMER0);

  build_map2(&mmers, mmer0_map, rlmap, mcmap, 
             mychunk, total_chunk, 
             mc_lower, mc_upper);

  int fd;
  struct stat sb;
  uint8_t *seq_p;

  fd = open(seqdb_file_path, O_RDONLY);
  if (fd == -1) handle_error("open");

  if (fstat(fd, &sb) == -1) /* To obtain file size */
    handle_error("fstat");

  khash_t(OVLP_CANDIDATES) * ovlp_candidates = kh_init(OVLP_CANDIDATES);
  build_ovlp_candidates(&mmers, mmer0_map, rlmap, mcmap, 
                        mychunk, total_chunk, 
                        mc_lower, mc_upper, 
                        ovlp_candidates);
  for (khiter_t __i = kh_begin(mmer0_map); __i != kh_end(mmer0_map); ++__i) {
    if (!kh_exist(mmer0_map, __i)) continue;
    mmer1_map = kh_val(mmer0_map, __i);
    for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
      if (!kh_exist(mmer1_map, __j)) continue;
      mpv = kh_val(mmer1_map, __j);
      kv_destroy(*mpv);
    }
    kh_destroy(MMER1, mmer1_map);
  }
  kh_destroy(MMER0, mmer0_map);

  kh_destroy(MMC, mcmap);
  kv_destroy(mmers);

  seq_p = (uint8_t *)mmap((void *)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
  dump_candidates(ovlp_candidates, rlmap, ovlp_upper, bestn, seq_p, ovlp_file);
  munmap(seq_p, sb.st_size);

  for (khiter_t __i = kh_begin(ovlp_candidates); __i != kh_end(ovlp_candidates); ++__i) {
    if (!kh_exist(ovlp_candidates, __i)) continue;
    ovlp_candidate_v *v = kh_val(ovlp_candidates, __i);
    kv_destroy(*v);
  }
  kh_destroy(OVLP_CANDIDATES, ovlp_candidates);

  kh_destroy(RLEN, rlmap);
  fclose(ovlp_file);
  if (!seqdb_prefix) free(seqdb_prefix);
  if (!shimmer_prefix) free(shimmer_prefix);
  if (!ovlp_file_path) free(ovlp_file_path);
}
