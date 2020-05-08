#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "kalloc.h"
#include "khash.h"
#include "kvec.h"
#include "shimmer.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define OVERLAP 0
#define CONTAINS 1
#define CONTAINED 2

#define handle_error(msg) \
  do {                    \
    perror(msg);          \
    exit(EXIT_FAILURE);   \
  } while (0)

#define SWAP(a, b)  {a ^= b; b ^= a; a ^= b;}; 

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
  uint32_t a_rid;
  uint32_t b_rid;
  int32_t match_size;
  double err_est;
  uint8_t a_strand;
  int32_t a_bgn, a_end;
  uint32_t a_len;
  uint8_t b_strand;
  int32_t b_bgn, b_end; 
  uint32_t b_len;
  char ovlp_type[24];
} ovlp_rec_t;

typedef struct {
  size_t n, m;
  ovlp_rec_t *a;
} ovlp_rec_v_t;

typedef struct {
  ovlp_rec_t ovlp;
  uint32_t m_bgn, m_end;
  hpc_seq_t * seq;
} ovlp_hpc_seq_t;

typedef struct {
  size_t n, m;
  ovlp_hpc_seq_t *a;
} ovlp_hpc_seq_v_t; 

typedef struct {
  uint64_t kmer_uint64;
  uint32_t pos;
} marker_t;

typedef struct {
  size_t n, m;
  marker_t *a;
} marker_v_t;


KHASH_MAP_INIT_INT(OVLP, ovlp_rec_v_t *);
KHASH_MAP_INIT_INT64(KMERCOUNT, uint32_t);

uint8_t twobit_map_f[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

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

seq_t * get_seq_from_db(uint8_t *seq_p,
                       uint32_t rid,
                       khash_t(RLEN) * rlmap,
                       uint8_t direction) {
  uint8_t *seq0;
  seq_t *seq;
  seq0 = get_read_seq_mmap_ptr(seq_p, rid, rlmap);
  rl_t rl;
  khiter_t k;
  k = kh_get(RLEN, rlmap, rid);
  assert(k != kh_end(rlmap));
  rl = kh_val(rlmap, k);
  seq = allocate_seq(rl.len+1);
  decode_biseq(seq0, seq->s, rl.len, direction);
  seq->l = rl.len;
  return seq;
}

hpc_seq_t * hp_compress(seq_t * seq) {
  assert(seq->l > 0);
  hpc_seq_t * cseq;  //homopolymer compressed seq
  char c = seq->s[0];
  size_t j = 0;
  cseq = allocate_hpc_seq(seq->l+1);
  cseq->s[0] = c;
  cseq->p[0] = 0;
  for (size_t i = 1; i < seq->l; i++) {
    if (seq->s[i] != c) {
      c = seq->s[i];
      j++;
      cseq->s[j] = c;
    }
    cseq->p[i] = j;
  }
  cseq->l = j;
  return cseq;
}

#define KMERSIZE 28
#define WINDOWSIZE 28
void fill_kmer_count(hpc_seq_t * seq, khash_t(KMERCOUNT) *kmer_count) {
  uint64_t mhash=0;
  mm128_t mmer0;
  khiter_t k;
  int32_t absent;

  mm128_v mmers = {0, 0, 0};
  kv_init(mmers);
  mm_sketch(NULL, seq->s, seq->l, WINDOWSIZE, KMERSIZE, 0, 0, &mmers);
  size_t s = 0;
  for (;;) {
    if (s >= mmers.n) break;
    mmer0 = mmers.a[s];
    mhash = mmer0.x >> 8;
    k = kh_put(KMERCOUNT, kmer_count, mhash, &absent);
    if (absent) {
      kh_value(kmer_count, k) = 1; 
    } else {
      kh_value(kmer_count, k) += 1; 
    }
    s++;
  }
  kv_destroy(mmers);
}  

void build_ovlp_map(FILE *ovlp_file, khash_t(OVLP) * ovlp_map,
                  uint32_t total_chunk, uint32_t my_chunk) {
  khiter_t k;
  uint32_t nr;
  ovlp_rec_t ovlp;
  ovlp_rec_v_t * ovlp_rec_v;
  int32_t absent;
  while (1) {
    nr = fscanf(ovlp_file, 
                "%d %d %d %lf %hhu %d %d %u %hhu %d %d %u %s\n",
                &(ovlp.a_rid), &(ovlp.b_rid), 
                &(ovlp.match_size), &(ovlp.err_est), 
                &(ovlp.a_strand), 
                &(ovlp.a_bgn), &(ovlp.a_end), &(ovlp.a_len),
                &(ovlp.b_strand), 
                &(ovlp.b_bgn), &(ovlp.b_end), &(ovlp.b_len),
                (char *) &(ovlp.ovlp_type));
    if (nr == 0 || nr == EOF) {
      break;
    }
    if (ovlp.a_rid % total_chunk == my_chunk % total_chunk) {
      k = kh_put(OVLP, ovlp_map, ovlp.a_rid, &absent);
      if (absent) {
        ovlp_rec_v = kmalloc(NULL, sizeof(ovlp_rec_v_t));
        ovlp_rec_v->n = 0;
        ovlp_rec_v->m = 0;
        ovlp_rec_v->a = NULL;
        kh_value(ovlp_map, k) = ovlp_rec_v;
      } else {
        ovlp_rec_v = kh_value(ovlp_map, k);
      }
      kv_push(ovlp_rec_t, NULL, *ovlp_rec_v, ovlp);
    }

    if (ovlp.b_rid % total_chunk == my_chunk % total_chunk) {
      k = kh_put(OVLP, ovlp_map, ovlp.b_rid, &absent);
      if (absent) {
        ovlp_rec_v = kmalloc(NULL, sizeof(ovlp_rec_v_t));
        ovlp_rec_v->n = 0;
        ovlp_rec_v->m = 0;
        ovlp_rec_v->a = NULL;
        kh_value(ovlp_map, k) = ovlp_rec_v;
      } else {
        ovlp_rec_v = kh_value(ovlp_map, k);
      }
      kv_push(ovlp_rec_t, NULL, *ovlp_rec_v, ovlp);
    }
  }
}
typedef khash_t(OVLP) OVLP_t;
typedef khash_t(RLEN) RLEN_t;
typedef khash_t(KMERCOUNT) KMERCOUNT_t;
void get_all_kmer_count(uint32_t rid, hpc_seq_t * a_cseq,
                        ovlp_rec_v_t *ovlp_rec_v, OVLP_t * ovlp_map,
                        uint8_t *seq_p, RLEN_t * rlmap,
                        uint32_t * cov, KMERCOUNT_t * kmer_count,
                        ovlp_hpc_seq_v_t *ovlp_seqs) {
  for (size_t i = 0; i < ovlp_rec_v->n; i++) {
    ovlp_rec_t ovlp;
    ovlp = ovlp_rec_v->a[i];
    uint32_t a_rid, b_rid;
    int32_t a_bgn, a_end;
    int32_t b_bgn, b_end;
    uint32_t a_len, b_len;
    a_rid = ovlp.a_rid;
    b_rid = ovlp.b_rid;
    a_bgn = ovlp.a_bgn >= 0 ? ovlp.a_bgn
                            : 0;  // maybe we need to fix negative a_bgn
    b_bgn = ovlp.b_bgn >= 0 ? ovlp.b_bgn : 0;
    a_end = ovlp.a_end;
    b_end = ovlp.b_end;
    a_len = ovlp.a_len;
    b_len = ovlp.b_len;
    if (ovlp.b_rid == rid) {
      SWAP(a_rid, b_rid);
      SWAP(a_bgn, b_bgn);
      SWAP(a_end, b_end);
      SWAP(a_len, b_len);
    }
    seq_t *b_seq;
    seq_t b_subseq;
    hpc_seq_t *b_csubseq;
    if (ovlp.b_strand == 0) {
      b_seq = get_seq_from_db(seq_p, b_rid, rlmap, ORIGINAL);
    } else {
      b_seq = get_seq_from_db(seq_p, b_rid, rlmap, REVERSED);
    }
    if (ovlp.b_strand == 1) {
      SWAP(b_bgn, b_end);
      b_bgn = b_len - b_bgn;
      b_end = b_len - b_end;
    }
    char *s1 = calloc(b_end - b_bgn + 1, sizeof(char));

    strncpy(s1, b_seq->s + b_bgn, b_end - b_bgn);
    b_subseq.l = b_end - b_bgn;
    b_subseq.m = b_subseq.l;
    b_subseq.s = s1;
    b_csubseq = hp_compress(&b_subseq);
    if (b_csubseq->l < 48) {
      free(s1);
      free_hpc_seq(b_csubseq);
      free_seq(b_seq);
      continue;
    };
    fill_kmer_count(b_csubseq, kmer_count);

    uint32_t m_bgn, m_end;
    m_bgn = a_cseq->p[a_bgn];
    m_end = a_cseq->p[a_end - 1];

    ovlp_hpc_seq_t ovlp_seq;
    ovlp_seq.m_bgn = m_bgn;
    ovlp_seq.m_end = m_end;
    ovlp_seq.ovlp = ovlp;
    ovlp_seq.seq = b_csubseq;
    kv_push(ovlp_hpc_seq_t, NULL, *ovlp_seqs, ovlp_seq);

    for (uint32_t i = m_bgn; i < m_end; i++) {
      cov[i] += 1;
    }

    free_seq(b_seq);
    free(s1);
  }
}

int main(int argc, char *argv[]) {
  char *seqdb_prefix = NULL;
  char *ovlp_file_path = NULL;

  FILE *ovlp_file;
  int fd;
  struct stat sb;
  //char ovlp_file_path[8192];
  char seq_idx_file_path[8192];
  char seqdb_file_path[8291];
  uint32_t total_chunk = 1, my_chunk = 1;
  int c;
  opterr = 0;
  while ((c = getopt(argc, argv, "p:l:t:c:")) != -1) {
    switch (c) {
      case 'p':
        seqdb_prefix = optarg;
        break;
      case 'l':
        ovlp_file_path = optarg;
        break;
      case 't':
        total_chunk = atoi(optarg);
        break;
      case 'c':
        my_chunk = atoi(optarg);
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
                  "Option -%c not specified, using 'preads.ovl' as the "
                  "overlap file\n",
                  optopt);
        }
        return 1;
      default:
        abort();
    }
  }

  assert(total_chunk > 0);
  assert(my_chunk > 0 && my_chunk <= total_chunk);
 
   if (seqdb_prefix == NULL) {
    seqdb_prefix = (char *)calloc(8192, 1);
    snprintf(seqdb_prefix, 8191, "seq_dataset");
  }

  if (ovlp_file_path == NULL) {
    ovlp_file_path = (char *)calloc(8192, 1);
    snprintf(ovlp_file_path, 8191, "preads.ovl");
  }
  
  int written;
  written = snprintf(seq_idx_file_path, sizeof(seq_idx_file_path), "%s.idx",
                     seqdb_prefix);
  assert(written < sizeof(seq_idx_file_path));
  fprintf(stderr, "using index file: %s\n", seq_idx_file_path);
  
  written = snprintf(seqdb_file_path, sizeof(seqdb_file_path), "%s.seqdb",
                     seqdb_prefix);
  assert(written < sizeof(seqdb_file_path));
  fprintf(stderr, "using seqdb file: %s\n", seqdb_file_path);
  fprintf(stderr, "using overlap file: %s\n", ovlp_file_path);
  
  fd = open(seqdb_file_path, O_RDONLY);
  if (fd == -1) handle_error("open");

  if (fstat(fd, &sb) == -1) /* To obtain file size */
    handle_error("fstat");

  khash_t(RLEN) * rlmap;
  rlmap = get_read_length_map(seq_idx_file_path);

  khash_t(OVLP) * ovlp_map = kh_init(OVLP);
  ovlp_rec_v_t * ovlp_rec_v;
  ovlp_file = fopen(ovlp_file_path, "r");
  build_ovlp_map(ovlp_file, ovlp_map, total_chunk, my_chunk);

  uint8_t *seq_p;
  seq_p = (uint8_t *)mmap((void *)0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
  // rlmap = get_read_length_map(seq_idx_file_path);
  seq_t *a_seq;
  hpc_seq_t *a_cseq;
  for (khiter_t __i = kh_begin(ovlp_map); __i != kh_end(ovlp_map); ++__i) {
    if (!kh_exist(ovlp_map, __i)) continue;

    ovlp_rec_v = kh_val(ovlp_map, __i);

    uint32_t rid;
    rid = kh_key(ovlp_map, __i);
    /*
    // useful for debugging
    if (rid != 1199) {
      continue;
    } 
    */
    khash_t(KMERCOUNT) *kmer_count = kh_init(KMERCOUNT);

    a_seq = get_seq_from_db(seq_p, rid, rlmap, ORIGINAL);
    a_cseq = hp_compress(a_seq);
    fill_kmer_count(a_cseq, kmer_count);

    ovlp_hpc_seq_v_t ovlp_seqs;
    kv_init(ovlp_seqs);

    uint32_t * cov;
    cov = calloc(a_cseq->l, sizeof(uint32_t));
    get_all_kmer_count(rid, a_cseq,
                       ovlp_rec_v, ovlp_map, seq_p, rlmap,
                       cov, kmer_count, &ovlp_seqs);

    uint32_t count = 0;
    marker_v_t markers;
    marker_t marker;
    uint32_t low_cov = 0;
    size_t p;
    uint64_t mhash=0;
    mm128_t mmer0;
    khiter_t k;

    kv_init(markers);

    for (p = 0; p < a_cseq->l; p++) {
      if (cov[p] < 2) {
        low_cov = 1;
        break;
      }
    }


    mm128_v mmers = {0, 0, 0};
    kv_init(mmers);
    mm_sketch(NULL, a_cseq->s, a_cseq->l, WINDOWSIZE, KMERSIZE, 0, 0, &mmers);
    size_t s = 0;
    for (;;) {
      if (s >= mmers.n) break;
      mmer0 = mmers.a[s];
      mhash = mmer0.x >> 8;
      p = (mmer0.y & 0xFFFFFFFF) >> 1;
      k = kh_get(KMERCOUNT, kmer_count, mhash);
      assert(k != kh_end(kmer_count));
      count = kh_val(kmer_count, k);

      if ((count > 4 && (double) count / (double) (cov[p]+1) > 0.1) && 
          (cov[p] - count > 3 && (double) count / (double) (cov[p]+1) < 0.9)) {
        marker.kmer_uint64 = mhash;
        marker.pos = p - KMERSIZE + 1;
        kv_push(marker_t, NULL, markers, marker);
      } 
      s++;
    }
    kv_destroy(mmers);

    for (size_t i=0; i < ovlp_seqs.n; i++) {
      ovlp_rec_t ovlp = ovlp_seqs.a[i].ovlp;
      if (rid != ovlp.a_rid) {
        free_hpc_seq(ovlp_seqs.a[i].seq);
        continue;
      }
      khash_t(KMERCOUNT) *seq_kmer_count = kh_init(KMERCOUNT);
      fill_kmer_count(ovlp_seqs.a[i].seq, seq_kmer_count); 
      uint32_t match_n = 0;
      uint32_t test_n = 0;
      uint32_t m_bgn, m_end;
      m_bgn = ovlp_seqs.a[i].m_bgn;
      m_end = ovlp_seqs.a[i].m_end;
      assert( m_bgn < m_end );
      for (size_t j=0; j < markers.n; j++) {
        marker_t m;
        khiter_t k;
        m = markers.a[j];
        if (m.pos >= m_bgn && m.pos < m_end - KMERSIZE + 1) {
          test_n += 1;
          k = kh_get(KMERCOUNT, seq_kmer_count, m.kmer_uint64);
          if (k != kh_end(seq_kmer_count)) {
            match_n += 1;
          }
        }
      }
      printf("%09d %09d %d %0.2f %hhu %d %d %u %hhu %d %d %u %s %d %d %d %d\n",
             ovlp.a_rid, ovlp.b_rid, ovlp.match_size, ovlp.err_est,
             ovlp.a_strand, ovlp.a_bgn, ovlp.a_end, ovlp.a_len, ovlp.b_strand,
             ovlp.b_bgn, ovlp.b_end, ovlp.b_len, (char *)ovlp.ovlp_type, test_n,
             match_n, test_n - match_n, low_cov);
      /*
      // useful for debugging
      if (ovlp.a_rid == 1199 && ovlp.b_rid == 468) {
        printf("seq: %s\n", ovlp_seqs.a[i].seq->s);
        printf("XX: %d %d\n", a_seq->l, a_cseq->l);
        printf(" Y: %d %d\n",m_bgn, m_end);
        for (size_t j=0; j < markers.n; j++) {
          marker_t m;
          khiter_t k;
          m = markers.a[j];
          if (m.pos >= m_bgn && m.pos < m_end ) {
            printf(" X:%d %lu\n", m.pos, m.kmer_uint64);
            test_n += 1;
            k = kh_get(KMERCOUNT, seq_kmer_count, m.kmer_uint64);
            if (k != kh_end(seq_kmer_count)) {
              printf("X2:%d %lu\n", m.pos, m.kmer_uint64);
              match_n += 1;
            }
          }
        }
      } 
      */
      kh_destroy(KMERCOUNT, seq_kmer_count);
      free_hpc_seq(ovlp_seqs.a[i].seq);
    }
    kv_destroy(ovlp_seqs);
    kv_destroy(markers);
    free_seq(a_seq);
    free_hpc_seq(a_cseq);
    free(cov);
    kh_destroy(KMERCOUNT, kmer_count);
  }

  for (khiter_t __i = kh_begin(ovlp_map); __i != kh_end(ovlp_map); ++__i) {
    if (!kh_exist(ovlp_map, __i)) continue;
    ovlp_rec_v = kh_val(ovlp_map, __i);
    kv_destroy(*ovlp_rec_v);
    free(ovlp_rec_v);
  }

  kh_destroy(OVLP, ovlp_map);
}
