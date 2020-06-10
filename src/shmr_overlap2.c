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
#define MMER_COUNT_UPPER_BOUND 240
#ifndef ORIGINAL
#define ORIGINAL 0
#endif
#ifndef REVERSED
#define REVERSED 1
#endif
#define READ_END_FUZZINESS 48
#define LOCAL_OVERLAP_UPPERBOUND 120
#define BESTN 4
#define OVERLAP 0
#define CONTAINS 1
#define CONTAINED 2
#define ALNBANDSIZE 100

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

mp128_v * match_shimmer_pair(khash_t(MMER0) * mmer0_map, mm128_t mmer0,mm128_t mmer1) {
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

KHASH_MAP_INIT_INT64(RPAIR, uint8_t);

bool check_rid_pair(khash_t(RPAIR) *rid_pairs, 
                    uint32_t rid0,
                    uint32_t rid1) {
  uint64_t ridp;
  khiter_t k;
  ridp = (((uint64_t)rid0) << 32) | ((uint64_t)rid1);
  k = kh_get(RPAIR, rid_pairs, ridp);
  return (k != kh_end(rid_pairs));
}

void set_rid_pair(khash_t(RPAIR) *rid_pairs,
                  uint32_t rid0,
                  uint32_t rid1) {
  uint64_t ridp;
  int32_t absent;
  khiter_t k;
  ridp = (((uint64_t)rid0) << 32) | ((uint64_t)rid1);
  k = kh_put(RPAIR, rid_pairs, ridp, &absent);
  kh_val(rid_pairs, k) = 1;
}

ovlp_match_t * match_seqs(uint8_t * seq0, uint8_t * seq1, 
                          uint32_t pos0, uint32_t pos1, 
                          uint32_t rlen0, uint32_t rlen1,
                          uint32_t strand0, uint32_t strand1) {
      uint32_t slen0;
      uint32_t slen1;
      ovlp_match_t *match;
      uint32_t align_bandwidth = 100;
      if (pos0 >= pos1) {
        slen0 = rlen0 - (pos0 - pos1);
        slen1 = rlen1;
        match = ovlp_match(seq0 + pos0 - pos1, slen0, strand0, 
                           seq1, slen1, strand1, 
                           align_bandwidth);
      } else {
        slen0 = rlen0;
        slen1 = rlen1 - (pos1 - pos0);
        match = ovlp_match(seq0, slen0, strand0, 
                           seq1 + pos1 - pos0, slen1, strand1, 
                           align_bandwidth);
      }
      return match;
}

void build_ovlp(mm128_v *mmers, 
               khash_t(MMER0) * mmer0_map,
               khash_t(RLEN) * rlmap, 
               khash_t(MMC) * mcmap, 
               uint32_t chunk,
               uint32_t total_chunk, 
               uint32_t lowerbound, 
               uint32_t upperbound,
               uint8_t * seq_p) {
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
  khash_t(RPAIR) *rid_pairs = kh_init(RPAIR);
  uint8_t *seq0 = NULL;
  uint8_t *seq1 = NULL;
  ovlp_match_t *match;

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
       mmer0 = mmer1;
       continue;
    }

      // don't use two minimers that are too close to each other
    if (((mmer1.y >> 1) & 0xFFFFFFF) - ((mmer0.y >> 1) & 0xFFFFFFF) < 100) {
      mmer0 = mmer1;
      continue;
    }

    mpv = match_shimmer_pair(mmer0_map, mmer0, mmer1);
    if (mpv == NULL) {
      mmer0 = mmer1;
      continue; 
    }

    rid1 = (uint32_t)(mmer0.y >> 32);
    pos1 = (uint32_t)((mmer0.y & 0xFFFFFFFF) >> 1) + 1;
    strand1 = ORIGINAL;
    k = kh_get(RLEN, rlmap, rid1);
    assert(k != kh_end(rlmap));
    rlen1 = kh_val(rlmap, k).len;
    for (size_t __k0 = 0; __k0 <  (mpv->n) ; __k0++) { 
      uint64_t y0;
      y0 = mpv->a[__k0].y0;
      rid0 = (uint32_t)(y0 >> 32);
      if (rid0 == rid1) {
          continue;
      }
      if (check_rid_pair(rid_pairs, rid0, rid1)) continue;

      k = kh_get(RLEN, rlmap, rid0);
      assert(k != kh_end(rlmap));
      rlen0 = kh_val(rlmap, k).len;

      pos0 = (uint32_t)((y0 & 0xFFFFFFFF) >> 1) + 1;
      strand0 = mpv->a[__k0].direction;

      seq0 = get_read_seq_mmap_ptr(seq_p, rid0, rlmap);
      seq1 = get_read_seq_mmap_ptr(seq_p, rid1, rlmap);
      match = match_seqs(seq0, seq1, 
                         pos0, pos1, 
                         rlen0, rlen1, 
                         strand0, strand1);
      
      
      seq_coor_t q_bgn, q_end, t_bgn, t_end;
      q_bgn = match->q_bgn;
      q_end = match->q_end;
      t_bgn = match->t_bgn;
      t_end = match->t_end;
      printf("%d %d %d %d %d %d %d %d %d %d\n", 
              rid0, pos0, strand0, 
              rid1, pos1, strand1,
              q_bgn, q_end,
              t_bgn, t_end);  
      free_ovlp_match(match);
      set_rid_pair(rid_pairs, rid0, rid1);
    }

    // reverse
    mpv = match_shimmer_pair(mmer0_map, mmer1, mmer0);
    if (mpv == NULL) {
      mmer0 = mmer1;
      continue; 
    }

    rid1 = (uint32_t)(mmer1.y >> 32);
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
      if (check_rid_pair(rid_pairs, rid0, rid1)) continue;

      k = kh_get(RLEN, rlmap, rid0);
      assert(k != kh_end(rlmap));
      rlen0 = kh_val(rlmap, k).len;

      pos0 = (uint32_t)((y0 & 0xFFFFFFFF) >> 1) + 1;
      strand0 = mpv->a[__k0].direction;

      seq0 = get_read_seq_mmap_ptr(seq_p, rid0, rlmap);
      seq1 = get_read_seq_mmap_ptr(seq_p, rid1, rlmap);
      match = match_seqs(seq0, seq1, 
                         pos0, pos1, 
                         rlen0, rlen1, 
                         strand0, strand1);

      seq_coor_t q_bgn, q_end, t_bgn, t_end;
      q_bgn = match->q_bgn;
      q_end = match->q_end;
      t_bgn = match->t_bgn;
      t_end = match->t_end;
      printf("%d %d %d %d %d %d %d %d %d %d\n", 
              rid0, pos0, strand0, 
              rid1, pos1, strand1,
              q_bgn, q_end,
              t_bgn, t_end);  
      set_rid_pair(rid_pairs, rid0, rid1);
      free_ovlp_match(match);
    }
    mmer0 = mmer1;
  }
  kh_destroy(RPAIR, rid_pairs);
}


int mp128_comp(const void *a, const void *b) {
  mp128_t *a0 = (mp128_t *)a;
  mp128_t *b0 = (mp128_t *)b;
  return ((a0->y0 & 0xFFFFFFFF) >> 1) < ((b0->y0 & 0xFFFFFFFF) >> 1);
}

void shimmer_to_overlap(mp128_v *mpv, khash_t(RLEN) * rlmap,
                        khash_t(RPAIR) * rid_pairs, uint8_t bestn,
                        uint32_t align_bandwidth, uint8_t *seq_p,
                        FILE *output) {
  uint64_t ridp;
  uint64_t y0;
  uint32_t rid0, pos0, rlen0, strand0;
  uint32_t rid1, pos1, rlen1, strand1;
  uint32_t right_ext = 0;
  khiter_t k;
  uint8_t *seq0 = NULL;
  uint8_t *seq1 = NULL;
  int32_t absent;
  uint8_t *contained;

  contained = (uint8_t *)calloc(
      mpv->n, sizeof(uint8_t));  // use calloc to set element to zero

  // clock_t time_begin = clock();
  // clock_t time_end;

  for (size_t __k0 = (mpv->n) - 1; __k0 > 0;
       __k0--) {  // note: k0 is an unsigned type
    if (contained[__k0 - 1] == 1) continue;
    y0 = mpv->a[__k0 - 1].y0;
    rid0 = (uint32_t)(y0 >> 32);
    pos0 = (uint32_t)((y0 & 0xFFFFFFFF) >> 1) + 1;
    k = kh_get(RLEN, rlmap, rid0);
    assert(k != kh_end(rlmap));
    rlen0 = kh_val(rlmap, k).len;
    strand0 = mpv->a[__k0 - 1].direction;
    seq0 = get_read_seq_mmap_ptr(seq_p, rid0, rlmap);

    if (right_ext == 0) {
      right_ext = rlen0;
    }

    size_t overlap_count = 0;
    for (size_t __k1 = 1; (__k0 + __k1 - 1 < mpv->n) && (overlap_count < bestn);
         __k1++) {
      if (contained[__k0 + __k1 - 1] == 1) continue;

      y0 = mpv->a[__k0 + __k1 - 1].y0;
      rid1 = (uint32_t)(y0 >> 32);

      if (rid0 == rid1) continue;
      // time_end = clock();
      // printf("X0: %lu %lu %lu %lu %09u %09u %lu\n", mpv->n, __k0, __k1,
      // overlap_count, rid0, rid1, time_end-time_begin); time_begin = time_end;
      ridp = rid0 < rid1 ? (((uint64_t)rid0) << 32) | ((uint64_t)rid1)
                         : (((uint64_t)rid1) << 32) | ((uint64_t)rid0);
      k = kh_get(RPAIR, rid_pairs, ridp);
      if (k != kh_end(rid_pairs)) {
        if (kh_val(rid_pairs, k) == OVERLAP) overlap_count += 1;
        continue;
      }
      pos1 = (uint32_t)((y0 & 0xFFFFFFFF) >> 1) + 1;
      assert(pos0 - pos1 >= 0);
      k = kh_get(RLEN, rlmap, rid1);
      assert(k != kh_end(rlmap));
      rlen1 = kh_val(rlmap, k).len;
      strand1 = mpv->a[__k0 + __k1 - 1].direction;
      seq1 = get_read_seq_mmap_ptr(seq_p, rid1, rlmap);

      // printf("X1: %lu %lu %lu %lu %09u %09u\n", mpv->n, __k0, __k1,
      // overlap_count, rid0, rid1);

      // printf("%09d %s\n%09d %s\n",rid0, seq0+pos0-pos1,rid1, seq1);
      uint32_t slen0 = rlen0 - pos0 + pos1;
      uint32_t slen1 = rlen1;
      ovlp_match_t *match;
      match = ovlp_match(seq0 + pos0 - pos1, slen0, strand0, seq1, slen1,
                         strand1, align_bandwidth);
      seq_coor_t q_bgn, q_end, t_bgn, t_end;
      q_bgn = match->q_bgn;
      q_end = match->q_end;
      t_bgn = match->t_bgn;
      t_end = match->t_end;
      // printf("X2: %u %u %d %d %d %d\n", pos0, pos1, q_bgn, q_end, t_bgn,
      // t_end);

      if ((q_bgn < READ_END_FUZZINESS && t_bgn < READ_END_FUZZINESS &&
           (abs((int64_t) slen0 - (int64_t) q_end) < READ_END_FUZZINESS ||
            abs((int64_t) slen1 - (int64_t) t_end) < READ_END_FUZZINESS)) &&
          q_end > 500 && t_end > 500) {
        // printf("X3: %u %u %d %d %d %d\n", pos0, pos1, q_bgn, q_end, t_bgn,
        // t_end); printf("%s\n%s\n", seq0+pos0-pos1, seq1);

        uint8_t ovlp_type;
        if (abs((int64_t) rlen0 - ( (int64_t) q_end - (int64_t) q_bgn)) < READ_END_FUZZINESS * 2 ||
            abs((int64_t) rlen1 - ( (int64_t) t_end - (int64_t) t_bgn)) < READ_END_FUZZINESS * 2) {
          if (rlen0 >= rlen1) {
            k = kh_put(RPAIR, rid_pairs, ridp, &absent);
            kh_val(rid_pairs, k) = CONTAINS;
            ovlp_type = CONTAINS;
            contained[__k0 + __k1 - 1] = 1;
          } else {
            k = kh_put(RPAIR, rid_pairs, ridp, &absent);
            kh_val(rid_pairs, k) = CONTAINED;
            ovlp_type = CONTAINED;
            contained[__k0 - 1] = 1;
          }
        } else {
          overlap_count++;
          k = kh_put(RPAIR, rid_pairs, ridp, &absent);
          kh_val(rid_pairs, k) = OVERLAP;
          ovlp_type = OVERLAP;
        }
        assert(absent == 1);

        ovlp_t ovlp;
        ovlp.y0 = mpv->a[__k0 - 1].y0;
        ovlp.y1 = mpv->a[__k0 + __k1 - 1].y0;
        ovlp.rl0 = rlen0;
        ovlp.rl1 = rlen1;
        ovlp.strand0 = strand0;
        ovlp.strand1 = strand1;
        ovlp.ovlp_type = ovlp_type;
        ovlp.match = *match;

        fwrite(&ovlp, sizeof(ovlp_t), 1, output);
      }
      free_ovlp_match(match);
      if (contained[__k0 - 1] == 1) break;
    }
  }
  free(contained);
}

void process_overlaps(char *seqdb_file_path, khash_t(MMER0) * mmer0_map,
                      khash_t(RLEN) * rlmap, khash_t(MMC) * mcmap,
                      uint8_t bestn, uint32_t ovlp_upper,
                      uint32_t align_bandwidth, FILE *ovlp_file) {
  int fd;
  struct stat sb;
  uint8_t *seq_p;
  mp128_v *mpv;
  // uint64_t mhash0, mhash1;

  khash_t(MMER1) * mmer1_map;

  fd = open(seqdb_file_path, O_RDONLY);
  if (fd == -1) handle_error("open");

  if (fstat(fd, &sb) == -1) /* To obtain file size */
    handle_error("fstat");

  seq_p = (uint8_t *)mmap((void *)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

  khash_t(RPAIR) *rid_pairs = kh_init(RPAIR);
  uint32_t iter_count = 0;
  clock_t begin = clock();
  clock_t end;
  for (khiter_t __i = kh_begin(mmer0_map); __i != kh_end(mmer0_map); ++__i) {
    if (!kh_exist(mmer0_map, __i)) continue;
    // mhash0 = kh_key(mmer0_map, __i);
    // mhash0 >>= 8;
    mmer1_map = kh_val(mmer0_map, __i);
    for (khiter_t __j = kh_begin(mmer1_map); __j != kh_end(mmer1_map); ++__j) {
      if (!kh_exist(mmer1_map, __j)) continue;
      // mhash1 = kh_key(mmer1_map, __j);
      // mhash1 >>= 8;
      mpv = kh_val(mmer1_map, __j);
      if (mpv->n <= 2 || mpv->n > ovlp_upper) continue;
      qsort(mpv->a, mpv->n, sizeof(mp128_t), mp128_comp);
      shimmer_to_overlap(mpv, rlmap, rid_pairs, bestn, align_bandwidth, seq_p,
                         ovlp_file);
      iter_count++;
      if (iter_count % 10000 == 0) {
        end = clock();
        fprintf(stderr, "iter_count: %d, %f s for 10000 candidates\n",
                iter_count, (double)(end - begin) / CLOCKS_PER_SEC);
        begin = end;
      }
    }
  }
  kh_destroy(RPAIR, rid_pairs);
  munmap(seq_p, sb.st_size);
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
  uint8_t bestn = BESTN;
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
  for (int i = 0; i < p.we_wordc; i++) {
    fprintf(stderr, "using shimmer count file: %s\n", mmc_fns[i]);
    mmc = read_mm_count(mmc_fns[i]);
    aggregate_mm_count(mcmap, &mmc);
    kv_destroy(mmc);
  }

  wordfree(&p);

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

  seq_p = (uint8_t *)mmap((void *)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

  build_ovlp(&mmers, mmer0_map, rlmap, mcmap, 
             mychunk, total_chunk, mc_lower,
             mc_upper,
             seq_p);
  munmap(seq_p, sb.st_size);
/*
  process_overlaps(seqdb_file_path, mmer0_map, rlmap, mcmap, bestn, ovlp_upper,
                   align_bandwidth, ovlp_file);
*/
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
  kh_destroy(RLEN, rlmap);
  kv_destroy(mmers);
  fclose(ovlp_file);
  if (!seqdb_prefix) free(seqdb_prefix);
  if (!shimmer_prefix) free(shimmer_prefix);
  if (!ovlp_file_path) free(ovlp_file_path);
}
