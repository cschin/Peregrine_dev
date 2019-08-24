from cffi import FFI
import os
import glob

basedir = os.environ["peregrine_base"]

ffibuilder = FFI()

ffibuilder.cdef("""
typedef struct {
    uint32_t max:31, zdropped:1;
    int max_q, max_t;      // max extension coordinate
    int mqe, mqe_t;        // max score when reaching the end of query
    int mte, mte_q;        // max score when reaching the end of target
    int score;             // max score reaching both ends; may be KSW_NEG_INF
    int m_cigar, n_cigar;
    int reach_end;
    uint32_t *cigar;
} ksw_extz_t;

void ksw_extz(void *km,
              int qlen, const uint8_t *query,
              int tlen, const uint8_t *target,
              int8_t m,
              const int8_t *mat,
              int8_t q, int8_t e, int w,
              int zdrop,
              int flag,
              ksw_extz_t *ez);

void ksw_extz2_sse(void *km,
                int qlen, const uint8_t *query,
                int tlen, const uint8_t *target,
                int8_t m,
                const int8_t *mat,
                int8_t q, int8_t e, int w,
                int zdrop,
                int end_bonus,
                int flag,
                ksw_extz_t *ez);

void ksw_extd(void *km,
                int qlen, const uint8_t *query,
                int tlen, const uint8_t *target,
                int8_t m, const int8_t *mat,
                int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2,
                int w, int zdrop,
                int flag,
                ksw_extz_t *ez);

void ksw_extd2_sse(void *km,
                   int qlen, const uint8_t *query,
                   int tlen, const uint8_t *target,
                   int8_t m, const int8_t *mat,
                   int8_t gapo, int8_t gape,
                   int8_t gapo2, int8_t gape2,
                   int w,
                   int zdrop,
                   int end_bonus,
                   int flag,
                   ksw_extz_t *ez);

void ksw_exts2_sse(void *km,
                   int qlen, const uint8_t *query,
                   int tlen, const uint8_t *target,
                   int8_t m, const int8_t *mat,
                   int8_t gapo, int8_t gape,
                   int8_t gapo2,
                   int8_t noncan,
                   int zdrop,
                   int flag,
                   ksw_extz_t *ez);

void ksw_extf2_sse(void *km,
                   int qlen, const uint8_t *query,
                   int tlen, const uint8_t *target,
                   int8_t mch, int8_t mis,
                   int8_t e, int w,
                   int xdrop, ksw_extz_t *ez);

ksw_extz_t * align(const char *tseq,
      const char *qseq,
      int sc_mch,
      int sc_mis,
      int gapo,
      int gape);
void free(void *ptr);
""")

srcs = glob.glob(f'{basedir}/ksw2/ksw2*.c') + [f'{basedir}/ksw2/align.c']

ffibuilder.set_source("peregrine._ksw4py",
               f"""
               #include "{basedir}/ksw2/ksw2.h"
               #include "{basedir}/ksw2/align.h"
               """, extra_compile_args = ["-march=native"],
               sources = srcs)   # library name, for the linker

if __name__ == "__main__":
    import sys
    ffibuilder.compile(verbose=True)
