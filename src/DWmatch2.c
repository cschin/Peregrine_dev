
/*
 *
 =====================================================================================
 *
 *       Filename:  DW_banded.c
 *
 *    Description:  A banded version for the O(ND) greedy sequence alignment
 algorithm
 *
 *        Version:  0.1
 *        Created:  07/20/2013 17:00:00
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jason Chin,
 *        Company:
 *
 *
 =====================================================================================

 #################################################################################$$
 # Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
 #
 # All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted (subject to the limitations in the
 # disclaimer below) provided that the following conditions are met:
 #
 #  * Redistributions of source code must retain the above copyright
 #  notice, this list of conditions and the following disclaimer.
 #
 #  * Redistributions in binary form must reproduce the above
 #  copyright notice, this list of conditions and the following
 #  disclaimer in the documentation and/or other materials provided
 #  with the distribution.
 #
 #  * Neither the name of Pacific Biosciences nor the names of its
 #  contributors may be used to endorse or promote products derived
 #  from this software without specific prior written permission.
 #
 # NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
 # GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
 # BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 # WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 # USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 # OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 # SUCH DAMAGE.
 #################################################################################$$
*/

#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "shimmer.h"

typedef struct {
  uint32_t l;  // length
  uint32_t m;  // allocated length
  char * s;  // sequence
} seq_t;

typedef struct {
  uint32_t l;  // length
  uint32_t m;  // allocated length
  char * s;  // sequence
  uint32_t * p;  // position map
} hpc_seq_t;

typedef struct {
  uint32_t x, y;
  int32_t pre_k;
} delta_point_t;

typedef struct {
  uint32_t n, m;
  delta_point_t *a;
} delta_point_v;


KHASH_MAP_INIT_INT64(BT, uint32_t);

seq_t * allocate_seq(uint32_t m) {
  seq_t * out;
  out = calloc(1, sizeof(seq_t));
  out->m = m;
  out->l = 0;
  out->s = calloc(m, sizeof(char));
  return out;
}

void free_seq(seq_t * seq) {
  free(seq->s);
  free(seq);
}

hpc_seq_t * allocate_hpc_seq(uint32_t m) {
  hpc_seq_t * out;
  out = calloc(1, sizeof(hpc_seq_t));
  out->m = m;
  out->l = 0;
  out->s = calloc(m, sizeof(char));
  out->p = calloc(m, sizeof(uint32_t));
  return out;
}

void free_hpc_seq(hpc_seq_t * cseq) {
  free(cseq->s);
  free(cseq->p);
  free(cseq);
}

hpc_seq_t * hp_compress(seq_t * seq) {
  assert(seq->l > 0);
  hpc_seq_t * cseq;  //homopolymer compressed seq
  char c = seq->s[0];
  size_t j = 0;
  cseq = allocate_hpc_seq(seq->l);
  cseq->s[0] = c;
  cseq->p[0] = 0;
  size_t i = 0;
  while (i < seq->l) {
    i ++;
    if (i < seq->l - 1 && i > 1 && seq->s[i-2] == seq->s[i] && seq->s[i-1] == seq->s[i+1]){
        i++;
        continue;
    }
    if (seq->s[i] != seq->s[i-1]) {
      c = seq->s[i];
      j++;
      cseq->s[j] = seq->s[i];
      cseq->p[j] = i;
    }
  }
  cseq->l = j;
  return cseq;
}

ovlp_match_t *ovlp_match2(uint8_t *query_seq, seq_coor_t q_len, uint8_t q_strand,
                         uint8_t *target_seq, seq_coor_t t_len,
                         uint8_t t_strand, seq_coor_t band_tolerance,
                         bool get_deltas,
                         bool use_hpc) {
  seq_coor_t *V;
  seq_coor_t *U;  // array of matched bases for each "k"
  seq_coor_t k_offset;
  seq_coor_t d;
  seq_coor_t k, k2;
  seq_coor_t best_m;  // the best "matches" for each d
  seq_coor_t min_k, new_min_k;
  seq_coor_t max_k, new_max_k;
  seq_coor_t x, y;
  seq_coor_t x1, y1;
  seq_coor_t max_d;
  seq_coor_t pre_k;
  seq_coor_t band_size;
  seq_t *qseq; // local queery seq
  seq_t *tseq; // local questy seq
  hpc_seq_t *qseq_hc;
  hpc_seq_t *tseq_hc;
  char *s0;
  char *s1;

  uint8_t q_shift = 0;
  uint8_t t_shift = 0;

  khash_t(BT) *bt_map;
  delta_point_v delta_pts;
  delta_point_t delta_pt;
  
  if (get_deltas == true) { 
   bt_map = kh_init(BT);
    kv_init(delta_pts);
  }

  qseq = allocate_seq(q_len);
  qseq->l = q_len;
  memcpy(qseq->s, query_seq, q_len);
  q_shift = q_strand == 0 ? 0 : 4;
  for (size_t i=0; i < q_len; i++) {
    qseq->s[i] = (qseq->s[i] >> q_shift) & 0x0F;
  }

  tseq = allocate_seq(t_len);
  tseq->l = t_len;
  memcpy(tseq->s, target_seq, t_len);
  t_shift = t_strand == 0 ? 0 : 4;
  for (size_t i=0; i < t_len; i++) {
    tseq->s[i] = (tseq->s[i] >> t_shift) & 0x0F;
  }
  
  if (use_hpc) {
    qseq_hc = hp_compress(qseq);
    tseq_hc = hp_compress(tseq);

    s0 = qseq_hc->s;
    s1 = tseq_hc->s;
    q_len = qseq_hc->l;
    t_len = tseq_hc->l;
  } else {
    s0 = qseq->s;
    s1 = tseq->s;
    q_len = qseq->l;
    t_len = qseq->l;
  }


  bool start = false;

  ovlp_match_t *rtn;
  bool matched = false;

  //printf("debug: %ld %ld\n", q_len, t_len);

  max_d = (int)(0.0015 * (q_len + t_len));

  band_size = band_tolerance * 2;

  V = calloc(max_d * 2 + 1, sizeof(seq_coor_t));
  U = calloc(max_d * 2 + 1, sizeof(seq_coor_t));

  k_offset = max_d;

  rtn = calloc(1, sizeof(ovlp_match_t));
  rtn->m_size = 0;
  rtn->q_bgn = 0;
  rtn->q_end = 0;
  rtn->t_bgn = 0;
  rtn->t_end = 0;
  rtn->q_m_end = 0;
  rtn->t_m_end = 0;
  rtn->reduce_deltas = NULL;
  uint32_t longest_match = 0;

  // printf("max_d: %lu, band_size: %lu\n", max_d, band_size);
  best_m = -1;
  min_k = 0;
  max_k = 0;
  for (d = 0; d < max_d; d++) {
    if (max_k - min_k > band_size) {
      break;
    }

    for (k = min_k; k <= max_k; k += 2) {
      if ((k == min_k) ||
          ((k != max_k) && (V[k - 1 + k_offset] < V[k + 1 + k_offset]))) {
        x = V[k + 1 + k_offset];
        pre_k = k +1;
      } else {
        x = V[k - 1 + k_offset] + 1;
        pre_k = k - 1;
      }
      y = x - k;

      if (get_deltas == true) {
        khiter_t _k;
        int32_t absent;
          // Note the k is signed int_32
        uint64_t key = ((uint64_t) d << 32) | (((uint64_t) k) & 0xFFFFFFFF);
        _k = kh_put(BT, bt_map, key, &absent);
        kh_val(bt_map, _k) = delta_pts.n;
        delta_pt.pre_k = pre_k;
        delta_pt.x = x;
        delta_pt.y = y;

        kv_push(delta_point_t, NULL, delta_pts, delta_pt);
      }

      x1 = x;
      y1 = y;


      while (x < q_len && y < t_len && s0[x] == s1[y]) {
        x++;
        y++;
      }

      if ((x - x1 > 8) && (start == false)) {
        if (use_hpc) {
          rtn->q_bgn = qseq_hc->p[x1];
          rtn->t_bgn = tseq_hc->p[y1];
        } else {
          rtn->q_bgn = x1;
          rtn->t_bgn = y1;
        }
        start = true;
      }

      if ((x - x1 > longest_match)) {
        longest_match = x - x1;
        if (use_hpc) {
          rtn->q_m_end = qseq_hc->p[x];
          rtn->t_m_end = tseq_hc->p[y];
        } else {
          rtn->q_m_end = x;
          rtn->t_m_end = y;
        }
      }

      V[k + k_offset] = x;
      U[k + k_offset] = x + y;

      if (x + y > best_m) {
        best_m = x + y;
      }

      if (x >= q_len || y >= t_len) {
        matched = true;
        break;
      }
    }

    // For banding
    new_min_k = max_k;
    new_max_k = min_k;

    for (k2 = min_k; k2 <= max_k; k2 += 2) {
      if (U[k2 + k_offset] >= best_m - band_tolerance) {
        if (k2 < new_min_k) {
          new_min_k = k2;
        }
        if (k2 > new_max_k) {
          new_max_k = k2;
        }
      }
    }

    max_k = new_max_k + 1;
    min_k = new_min_k - 1;

    if (matched == true) {
      //printf("match %d %d %d %d\n", x, y, qseq_hc->p[x], tseq_hc->p[y]);
      if (use_hpc) {
        rtn->q_end = qseq_hc->p[x];
        rtn->t_end = tseq_hc->p[y];
      } else {
        rtn->q_end = x;
        rtn->t_end = y;
      }
      rtn->dist = d;
      // we don't really generate the alingment path here, so we can only
      // estimate the alignment string size
      rtn->m_size =
          (rtn->q_end - rtn->q_bgn + rtn->t_end - rtn->t_bgn + 2 * d) / 2;

      if (get_deltas == true && d > 0) {
        //printf("== %d %d\n", d, k);
        rtn->reduce_deltas = calloc(d, sizeof(uint32_t)); 
        uint32_t dd;
        int32_t kk;
        dd = d;
        kk = k;
        while (dd > 0 || kk != 0) {
          khiter_t _k;
          // Note the kk is signed int_32
          uint64_t key = ((uint64_t) dd << 32) | (((uint64_t) kk) & 0xFFFFFFFF);
          _k = kh_get(BT, bt_map, key);
          assert (_k != kh_end(bt_map));
          uint32_t n = kh_val(bt_map, _k);
          delta_pt = delta_pts.a[n];
          if (use_hpc) {
            uint32_t xx = qseq_hc->p[delta_pt.x];
            uint32_t yy = qseq_hc->p[delta_pt.y];

            if (delta_pt.pre_k - kk > 0) {
              rtn->reduce_deltas[dd-1] = (xx << 1); 
            } else {
              rtn->reduce_deltas[dd-1] = (xx << 1) | 0x1; 
            }

          } else {
            if (delta_pt.pre_k - kk > 0) {
              rtn->reduce_deltas[dd-1] = (x << 1); 
            } else {
              rtn->reduce_deltas[dd-1] = (x << 1) | 0x1; 
            }
          }
          dd--;
          kk = delta_pt.pre_k;
        }
      }
      break;
    }
  }
  if (matched == false) {
    rtn->q_bgn = 0;
    rtn->t_bgn = 0;
    rtn->m_size = 0;
    rtn->dist = 1;
  }

  free_seq(qseq);
  free_seq(tseq);
  if (use_hpc) {
    free_hpc_seq(qseq_hc);
    free_hpc_seq(tseq_hc);
  }
  free(V);
  free(U);
  if (get_deltas == true) {
    kv_destroy(delta_pts);
    kh_destroy(BT, bt_map);
  }
  return rtn;
}

void free_ovlp_match(ovlp_match_t *match) { 
  if (match->reduce_deltas != NULL) {
    free(match->reduce_deltas);
  }
  free(match); 
}
