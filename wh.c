/*
 * Copyright (c) 2016--2018  Wu, Xingbo <wuxb45@gmail.com>
 *
 * All rights reserved. No warranty, explicit or implicit, provided.
 */
#define _GNU_SOURCE

// headers {{{
#include "wh.h"
#include <sched.h>
#include <assert.h>
#include <math.h>
#include <execinfo.h>
#include <signal.h>
#include <stdatomic.h>
#include <sys/socket.h>
#include <sys/sysinfo.h>
#include <netdb.h>
#include <byteswap.h>
#include <x86intrin.h>
// }}} headers

// ansi colors {{{
#define ANSI_FR  "\x1b[31m"

#define ANSI_X   "\x1b[0m"
// }}} ansi colors

// atomic {{{
/* C11 atomic types */
typedef atomic_bool             abool;

typedef atomic_uint_least8_t    au8;
typedef atomic_uint_least16_t   au16;
typedef atomic_uint_least32_t   au32;
typedef atomic_uint_least64_t   au64;

typedef atomic_int_least8_t     as8;
typedef atomic_int_least16_t    as16;
typedef atomic_int_least32_t    as32;
typedef atomic_int_least64_t    as64;
// }}} atomic

// locking {{{
  inline void
spinlock_init(spinlock * const lock)
{
  lock->var = 0u;
}

  inline void
spinlock_lock(spinlock * const lock)
{
  au16 * const pvar = (typeof(pvar))&(lock->var);
  do {
    if (0u == atomic_fetch_add(pvar, 1u))
      return;
    do {
      _mm_pause();
    } while (lock->var);
  } while (true);
}

  inline bool
spinlock_trylock(spinlock * const lock)
{
  au16 * const pvar = (typeof(pvar))&(lock->var);
  u16 z = 0u;
  return atomic_compare_exchange_weak(pvar, &z, 1u);
}

  inline void
spinlock_unlock(spinlock * const lock)
{
  lock->var = 0u;
}

  inline void
mutexlock_init(mutexlock * const lock)
{
  pthread_mutex_init(&(lock->lock), NULL);
}

  inline void
mutexlock_lock(mutexlock * const lock)
{
  do {
    const int r = pthread_mutex_lock(&(lock->lock));
    if (r == 0) return;
    else if (r != EAGAIN) exit(0);
  } while (true);
}

  inline bool
mutexlock_trylock(mutexlock * const lock)
{
  do {
    const int r = pthread_mutex_trylock(&(lock->lock));
    if (r == 0) return true;
    else if (r == EBUSY) return false;
    else if (r != EAGAIN) exit(0);
  } while (true);
}

  inline void
mutexlock_unlock(mutexlock * const lock)
{
  do {
    const int r = pthread_mutex_unlock(&(lock->lock));
    if (r == 0) return;
    else if ((r != EAGAIN)) exit(0);
  } while (true);
}

#define RWLOCK_WSHIFT ((15))
#define RWLOCK_WBIT ((1u << RWLOCK_WSHIFT))

  inline void
rwlock_init(rwlock * const lock)
{
  au16 * const pvar = (typeof(pvar))(&(lock->var));
  atomic_store(pvar, 0);
}

  inline bool
rwlock_trylock_read(rwlock * const lock)
{
  au16 * const pvar = (typeof(pvar))(&(lock->var));
  if ((atomic_fetch_add(pvar, 1u) >> RWLOCK_WSHIFT) == 0u) {
    return true;
  } else {
    atomic_fetch_sub(pvar, 1u);
    return false;
  }
}

// actually nr + 1
  inline bool
rwlock_trylock_read_nr(rwlock * const lock, u64 nr)
{
  if (rwlock_trylock_read(lock))
    return true;
  do {
    if (((lock->var >> RWLOCK_WSHIFT) == 0u) && rwlock_trylock_read(lock))
      return true;
    _mm_pause();
  } while (nr--);
  return false;
}

  inline void
rwlock_lock_read(rwlock * const lock)
{
  while (rwlock_trylock_read(lock) == false)
    while (lock->var >> RWLOCK_WSHIFT)
      _mm_pause();
}

  inline void
rwlock_unlock_read(rwlock * const lock)
{
  au16 * const pvar = (typeof(pvar))(&(lock->var));
  atomic_fetch_sub(pvar, 1u);
}

  inline bool
rwlock_trylock_write(rwlock * const lock)
{
  au16 * const pvar = (typeof(pvar))(&(lock->var));
  u16 v0 = *pvar;
  if (v0 == 0u) {
    if (atomic_compare_exchange_weak(pvar, &v0, RWLOCK_WBIT))
      return true;
  }
  return false;
}

// actually nr + 1
  inline bool
rwlock_trylock_write_nr(rwlock * const lock, u64 nr)
{
  do {
    if (rwlock_trylock_write(lock))
      return true;
    _mm_pause();
  } while (nr--);
  return false;
}

  inline void
rwlock_lock_write(rwlock * const lock)
{
  while (rwlock_trylock_write(lock) == false)
    while (lock->var)
      _mm_pause();
}

  inline void
rwlock_unlock_write(rwlock * const lock)
{
  au16 * const pvar = (typeof(pvar))(&(lock->var));
  atomic_fetch_sub(pvar, RWLOCK_WBIT);
}
#undef RWLOCK_WBIT
// }}} locking

// timing {{{
  inline u64
rdtsc(void)
{
  return _rdtsc();
}

  inline u64
time_nsec(void)
{
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
  return ts.tv_sec * UINT64_C(1000000000) + ts.tv_nsec;
}

  inline double
time_sec(void)
{
  const u64 nsec = time_nsec();
  return ((double)nsec) / 1000000000.0;
}

  inline u64
time_diff_nsec(const u64 last)
{
  return time_nsec() - last;
}

  inline double
time_diff_sec(const double last)
{
  return time_sec() - last;
}

  inline u64
timespec_diff(const struct timespec t0, const struct timespec t1)
{
  return UINT64_C(1000000000) * (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec);
}
// }}} timing

// debug {{{
  void
debug_break(void)
{
  usleep(100);
}

  void
debug_backtrace(void)
{
  void *array[100];
  const int size = backtrace(array, 100);
  dprintf(2, "Backtrace (%d):\n", size);
  // skip this call
  backtrace_symbols_fd(array + 1, size, 2);
}

  void
debug_wait_gdb(void)
{
  debug_backtrace();
  bool wait = true;
  volatile bool * const v = &wait;
  *v = true;

  time_t now;
  time(&now);
  struct tm nowtm;
  localtime_r(&now, &nowtm);
  char timestamp[64] = {};
  strftime(timestamp, 64, "%F %T %Z (%z)", &nowtm);

  char hostname[256] = {};
  gethostname(hostname, 256);
  char threadname[256];
  pthread_getname_np(pthread_self(), threadname, 256);

  const char * const pattern = "[Waiting GDB] %s %s @ %s\n"
    "    Attach me:   " ANSI_FR "sudo -Hi gdb -p %d" ANSI_X "\n";
  fprintf(stderr, pattern, timestamp, threadname, hostname, getpid());
  fflush(stderr);
  // to continue: gdb> set var *v = 0
  while (*v) {
    sleep(1);
  }
}

#ifndef NDEBUG
  inline void
debug_assert(const bool v)
{
  if (!v) debug_wait_gdb();
}
#endif

__attribute__((noreturn))
  inline void
debug_die(void)
{
  debug_wait_gdb();
  exit(0);
}

  static void
__signal_handler_wait_gdb(const int sig, siginfo_t * const info, void * const context)
{
  (void)info;
  (void)context;
  printf("[SIGNAL] %s\n", strsignal(sig));
  debug_wait_gdb();
}

__attribute__((constructor))
  static void
debug_catch_fatals(void)
{
  struct sigaction sa = {};
  sa.sa_sigaction = __signal_handler_wait_gdb;
  sigemptyset(&(sa.sa_mask));
  sa.sa_flags = SA_SIGINFO;
  const int fatals[] = {SIGSEGV, SIGFPE, SIGILL, SIGBUS, 0};
  for (int i = 0; fatals[i]; i++) {
    if (sigaction(fatals[i], &sa, NULL) == -1) {
      fprintf(stderr, "Failed to set signal handler for %s\n", strsignal(fatals[i]));
      fflush(stderr);
    }
  }
}
// }}} debug

// process/thread {{{
  u64
process_get_rss(void)
{
  u64 size;
  u64 rss;
  FILE * const fp = fopen("/proc/self/statm", "r");
  if (fp == NULL) return 0;
  if (fscanf(fp, "%" PRIu64 "%" PRIu64, &size, &rss) != 2) {
    fclose(fp);
    return 0;
  }
  fclose(fp);
  return rss * (u64)sysconf(_SC_PAGESIZE);
}

  static inline cpu_set_t *
__cpu_set_alloc(size_t * const size_out)
{
  const int ncpu = sysconf(_SC_NPROCESSORS_CONF);
  const size_t s1 = CPU_ALLOC_SIZE(ncpu);
  const size_t s2 = sizeof(cpu_set_t);
  const size_t size = s1 > s2 ? s1 : s2;
  *size_out = size;
  cpu_set_t * const set = malloc(size);
  return set;
}

  u64
process_affinity_core_count(void)
{
  size_t xsize = 0;
  cpu_set_t * const set = __cpu_set_alloc(&xsize);
  if (sched_getaffinity(0, xsize, set) != 0) {
    return sysconf(_SC_NPROCESSORS_CONF);
  }
  const int nr = CPU_COUNT_S(xsize, set);
  free(set);
  return (nr > 0) ? nr : sysconf(_SC_NPROCESSORS_CONF);
}

  u64
process_affinity_core_list(const u64 max, u64 * const cores)
{
  memset(cores, 0, max * sizeof(cores[0]));
  size_t xsize = 0;
  cpu_set_t * const set = __cpu_set_alloc(&xsize);
  if (sched_getaffinity(0, xsize, set) != 0) return 0;

  const u64 ncpu = sysconf(_SC_NPROCESSORS_CONF);
  const u64 nr_affinity = CPU_COUNT_S(xsize, set);
  const u64 nr = nr_affinity < max ? nr_affinity : max;
  u64 j = 0;
  for (u64 i = 0; i < ncpu; i++) {
    if (CPU_ISSET_S((int)i, xsize, set)) {
      cores[j++] = i;
    }
    if (j >= nr) break;
  }
  free(set);
  return j;
}

  u64
process_cpu_time_usec(void)
{
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  const u64 usr = (r.ru_utime.tv_sec * UINT64_C(1000000)) + r.ru_utime.tv_usec;
  const u64 sys = (r.ru_stime.tv_sec * UINT64_C(1000000)) + r.ru_stime.tv_usec;
  return usr + sys;
}

  void
thread_set_affinity(const u64 cpu)
{
  const u64 ncpu = sysconf(_SC_NPROCESSORS_CONF);
  size_t xsize = 0;
  cpu_set_t * const set = __cpu_set_alloc(&xsize);

  CPU_ZERO_S(xsize, set);
  CPU_SET_S(cpu % ncpu, xsize, set);
  sched_setaffinity(0, xsize, set);
  free(set);
}

  double
thread_fork_join_private(const u64 nr, void *(*func) (void *), void * const * const argv)
{
  if (nr == 0) return 0.0;
  const u64 ncpu = sysconf(_SC_NPROCESSORS_CONF);
  int cores[ncpu];
  size_t xsize = 0;
  cpu_set_t * const set = __cpu_set_alloc(&xsize);

  const bool force_all = (sched_getaffinity(0, xsize, set) != 0) ? true : false;
  u64 j = 0;
  for (u64 i = 0; i < ncpu; i++) {
    if (force_all || CPU_ISSET_S((int)i, xsize, set)) {
      cores[j++] = i;
    }
  }
  const u64 ncores = j;

  const u64 nr_threads = nr ? nr : ncores;
  const double t0 = time_sec();
  if (nr_threads == 1lu) { // no fork for one thread
    size_t xsize0 = 0;
    const pthread_t self = pthread_self();
    cpu_set_t * const set0 = __cpu_set_alloc(&xsize0);
    pthread_getaffinity_np(self, xsize0, set0);
    CPU_ZERO_S(xsize, set);
    CPU_SET_S(cores[0], xsize, set);
    pthread_setaffinity_np(self, xsize, set);
    func(argv[0]);
    pthread_setaffinity_np(self, xsize0, set0);
    free(set0);
  } else { // fork
    pthread_t tids[nr_threads];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    char thname[32];
    for (u64 i = 0; i < nr_threads; i++) {
      CPU_ZERO_S(xsize, set);
      CPU_SET_S(cores[i % ncores], xsize, set);
      pthread_attr_setaffinity_np(&attr, xsize, set);
      const int r = pthread_create(&(tids[i]), &attr, func, argv[i]);
      if (r != 0) {
        tids[i] = 0;
      } else {
        sprintf(thname, "fork_join_%"PRIu64, i);
        pthread_setname_np(tids[i], thname);
      }
    }
    pthread_attr_destroy(&attr);
    for (u64 i = 0; i < nr_threads; i++) {
      if (tids[i]) pthread_join(tids[i], NULL);
    }
  }
  const double dt = time_diff_sec(t0);
  free(set);
  return dt;
}

  inline double
thread_fork_join(const u64 nr, void *(*func) (void *), void * const arg)
{
  const u64 nthreads = nr ? nr : process_affinity_core_count();
  void * argv[nthreads];
  for (u64 i = 0; i < nthreads; i++) {
    argv[i] = arg;
  }
  return thread_fork_join_private(nthreads, func, argv);
}

  inline int
thread_create_at(const u64 cpu, pthread_t * const thread, void *(*start_routine) (void *), void * const arg)
{
  const u64 ncpu = sysconf(_SC_NPROCESSORS_CONF);
  const u64 cpu_id = cpu % ncpu;
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  size_t xsize = 0;
  cpu_set_t * const set = __cpu_set_alloc(&xsize);

  CPU_ZERO_S(xsize, set);
  CPU_SET_S(cpu_id, xsize, set);
  pthread_attr_setaffinity_np(&attr, xsize, set);
  const int r = pthread_create(thread, &attr, start_routine, arg);
  pthread_attr_destroy(&attr);
  free(set);
  return r;
}

  inline u64
thread_get_core(void)
{
  return (u64)sched_getcpu();
}
// }}} process/thread

// mm {{{
#define PGSZ ((UINT64_C(4096)))
  inline void *
xalloc(const u64 align, const u64 size)
{
  void * p;
  const int r = posix_memalign(&p, align, size);
  if (r == 0) return p;
  else return NULL;
}

// alloc cache-line aligned address
  inline void *
yalloc(const u64 size)
{
  void * p;
  const int r = posix_memalign(&p, 64, size);
  if (r == 0) return p;
  else return NULL;
}

  inline void
pages_unmap(void * const ptr, const size_t size)
{
#ifndef HEAPCHECKING
  munmap(ptr, size);
#else
  (void)size;
  free(ptr);
#endif
}

  static inline void *
__pages_alloc(const size_t size, const int flags)
{
  // vi /etc/security/limits.conf
  // * - memlock unlimited
  void * const p = mmap(NULL, size, PROT_READ | PROT_WRITE, flags, -1, 0);
  if (p == MAP_FAILED) {
    return NULL;
  }
  mlock(p, size); // ignore if cannot pin memory. see memlock in /etc/security/limits.conf
  return p;
}

  inline void *
pages_alloc_1gb(const size_t nr_1gb)
{
  const u64 sz = nr_1gb << 30;
#ifndef HEAPCHECKING
  return __pages_alloc(sz, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | (30 << MAP_HUGE_SHIFT));
#else
  void * const p = xalloc(UINT64_C(1) << 30, sz);
  if (p) memset(p, 0, sz);
  return p;
#endif
}

  inline void *
pages_alloc_2mb(const size_t nr_2mb)
{
  const u64 sz = nr_2mb << 21;
#ifndef HEAPCHECKING
  return __pages_alloc(sz, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | (21 << MAP_HUGE_SHIFT));
#else
  void * const p = xalloc(UINT64_C(1) << 21, sz);
  if (p) memset(p, 0, sz);
  return p;
#endif
}

  inline void *
pages_alloc_4kb(const size_t nr_4kb)
{
  const size_t sz = nr_4kb << 12;
#ifndef HEAPCHECKING
  return __pages_alloc(sz, MAP_PRIVATE | MAP_ANONYMOUS);
#else
  void * const p = xalloc(UINT64_C(1) << 12, sz);
  if (p) memset(p, 0, sz);
  return p;
#endif
}

  void *
pages_alloc_best(const size_t size, const bool try_1gb, u64 * const size_out)
{
  if (try_1gb) {
    const size_t nr_1gb = (size + ((UINT64_C(1) << 30) - UINT64_C(1))) >> 30;
    // 1gb super huge page: waste < 1/16 or 6.25%
    if (((nr_1gb << 30) - size) < (size >> 4)) {
      void * const p1 = pages_alloc_1gb(nr_1gb);
      if (p1) {
        *size_out = nr_1gb << 30;
        return p1;
      }
    }
  }

  // 2mb huge page: at least 1MB
  if (size >= (UINT64_C(1) << 20)) {
    const size_t nr_2mb = (size + ((UINT64_C(1) << 21) - UINT64_C(1))) >> 21;
    void * const p2 = pages_alloc_2mb(nr_2mb);
    if (p2) {
      *size_out = nr_2mb << 21;
      return p2;
    }
  }

  const size_t nr_4kb = (size + ((UINT64_C(1) << 12) - UINT64_C(1))) >> 12;
  void * const p3 = pages_alloc_4kb(nr_4kb);
  if (p3) {
    *size_out = nr_4kb << 12;
  }
  return p3;
}
// }}} mm

// slab {{{
struct slab_object {
  struct slab_object * next;
};

struct slab {
  spinlock lock;
  u64 obj_size;
  struct slab_object * obj_head;
  u64 blk_size; // size of each memory block
  u64 blk_nr; // number of blocks
  u64 nr_alloc; // number of objects in use
  struct slab_object * blk_head; // list of all blocks
};

  static void
slab_expand(struct slab * const slab)
{
  size_t blk_size;
  struct slab_object * const blk = pages_alloc_best(slab->blk_size, true, &blk_size);
  (void)blk_size;
  debug_assert(blk_size == slab->blk_size);
  blk->next = slab->blk_head;
  slab->blk_head = blk;

  const u64 nr_objs = (slab->blk_size - 64u) / slab->obj_size;
  for (u64 i = nr_objs; i; i--) {
    struct slab_object * const obj = (typeof(obj))(((u8 *)blk) + 64u + ((i - 1lu) * slab->obj_size));
    obj->next = slab->obj_head;
    slab->obj_head = obj;
  }
}

  struct slab *
slab_create(const u64 obj_size, const u64 blk_size)
{
  // obj must be 8-byte aligned
  // blk must be at least of page size and power of 2
  if ((obj_size % 8lu) || (blk_size < 4096lu) || (blk_size & (blk_size - 1lu)))
    return NULL;

  struct slab * const slab = (typeof(slab))malloc(sizeof(*slab));
  debug_assert(slab);
  spinlock_init(&(slab->lock));
  slab->obj_size = obj_size;
  slab->obj_head = NULL;
  slab->blk_size = blk_size;
  slab->blk_nr = 0lu;
  slab->blk_head = NULL;
  slab->nr_alloc = 0lu;
  slab_expand(slab);
  return slab;
}

  inline void *
slab_alloc_unsafe(struct slab * const slab)
{
  if (slab->obj_head == NULL) {
    slab_expand(slab);
  }
  debug_assert(slab->obj_head);
  struct slab_object * const obj = slab->obj_head;
  slab->obj_head = obj->next;
  slab->nr_alloc++;
  return (void *)obj;
}

  inline void *
slab_alloc(struct slab * const slab)
{
  spinlock_lock(&(slab->lock));
  void * const ptr = slab_alloc_unsafe(slab);
  spinlock_unlock(&(slab->lock));
  return ptr;
}

  inline void
slab_free_unsafe(struct slab * const slab, void * const ptr)
{
  struct slab_object * const obj = (typeof(obj))ptr;
  obj->next = slab->obj_head;
  slab->obj_head = obj;
  slab->nr_alloc--;
}

  inline void
slab_free(struct slab * const slab, void * const ptr)
{
  spinlock_lock(&(slab->lock));
  slab_free_unsafe(slab, ptr);
  spinlock_unlock(&(slab->lock));
}

  void
slab_destroy(struct slab * const slab)
{
  struct slab_object * iter = slab->blk_head;
  while (iter) {
    struct slab_object * const next = iter->next;
    pages_unmap(iter, slab->blk_size);
    iter = next;
  }
  free(slab);
}
// }}} slab

// cpucache {{{
  inline void
cpu_mfence(void)
{
  atomic_thread_fence(memory_order_seq_cst);
}

// compiler fence
  inline void
cpu_cfence(void)
{
  atomic_thread_fence(memory_order_acq_rel);
}
// }}} cpucache

// hash {{{
#define CRC32C_SEED ((0xDEEDBEEFu))

// for crc less than 3 bytes
  static inline u32
crc32c_inc_short(const u8 * buf, size_t nr, u32 crc)
{
  debug_assert(nr <= 3);
  if (nr == 0)
    return crc;

  // nr == 1
  crc = _mm_crc32_u8(crc, *buf);
  if ((--nr) == 0)
    return crc;

  buf++;
  // nr == 2
  crc = _mm_crc32_u8(crc, *buf);
  if ((--nr) == 0)
    return crc;

  buf++;
  // nr == 3
  return _mm_crc32_u8(crc, *buf);
}

  static inline u32
crc32c_inc(const u8 * buf, size_t nr, u32 crc)
{
#pragma nounroll
  while (nr >= sizeof(u64)) {
    crc = _mm_crc32_u64(crc, *((u64*)buf));
    nr -= sizeof(u64);
    buf += sizeof(u64);
  }
  if (nr >= sizeof(u32)) {
    crc = _mm_crc32_u32(crc, *((u32*)buf));
    nr -= sizeof(u32);
    buf += sizeof(u32);
  }
  return crc32c_inc_short(buf, nr, crc);
}

  inline u32
crc32c(const void * const ptr, size_t len)
{
  return crc32c_inc((const u8 *)ptr, len, CRC32C_SEED);
}

  static inline u64
crc32c_extend(const u32 lo)
{
  const u64 hi = (u64)(~lo);
  return (hi << 32) | ((u64)lo);
}
// }}} hash

// rcu {{{

// node {{{
// bits 63 -- 16 value (pointer)
// bits 15       valid
// bits 14 -- 0  count (refcount)
#define RCU_COUNT_MASK   ((UINT64_C(0x0000000000007fff)))
#define RCU_VALID_MASK   ((UINT64_C(0x0000000000008000)))
#define RCU_VALUE_MASK   ((UINT64_C(0xffffffffffff0000)))
#define RCU_VALUE_SHIFT  ((16))
struct rcu_node {
  au64 x[8];
};

  void
rcu_node_init(struct rcu_node * const node)
{
  atomic_store(&(node->x[0]), RCU_VALID_MASK); // valid null pointer
  atomic_store(&(node->x[1]), RCU_VALUE_MASK); // invalid non-null pointer
}

  struct rcu_node *
rcu_node_create(void)
{
  struct rcu_node * const node = (typeof(node))malloc(sizeof(*node));
  rcu_node_init(node);
  return node;
}

  void *
rcu_node_ref(struct rcu_node * const node)
{
  do {
    for (u64 i = 0; i < 2; i++) {
      u64 x = node->x[i];
      if (x & RCU_VALID_MASK) {
        if (atomic_compare_exchange_weak(&(node->x[i]), &x, x + 1))
          return (void *)(x >> RCU_VALUE_SHIFT);
      }
    }
  } while (true);
}

  void
rcu_node_unref(struct rcu_node * const node, void * const ptr)
{
  for (u64 i = 0; i < 2; i++) {
    const u64 x = node->x[i];
    if ((x >> RCU_VALUE_SHIFT) == ((u64)ptr)) {
      (void)atomic_fetch_sub(&(node->x[i]), 1lu);
      return;
    }
  }
  debug_die();
}

  void
rcu_node_update(struct rcu_node * const node, void * const ptr)
{
  const u64 xx = (((u64)ptr) << RCU_VALUE_SHIFT) | RCU_VALID_MASK;
  for (u64 i = 0; i < 2; i++) {
    if ((node->x[i] & RCU_VALID_MASK) == 0) {
      atomic_store(&(node->x[i]), xx);
      au64 * const v = &(node->x[1-i]);
      debug_assert(((*v) >> RCU_VALUE_SHIFT) != ((u64)ptr));
      (void)atomic_fetch_sub(v, RCU_VALID_MASK);
      while (atomic_load(v) & RCU_COUNT_MASK)
        _mm_pause();
      return;
    }
  }
  debug_die();
}
// }}} node

// multi-rcu {{{
struct rcu {
  u64 nr;
  u64 mask;
  u64 padding[6];
  struct rcu_node nodes[];
};

  void
rcu_init(struct rcu * const rcu, const u64 nr)
{
  rcu->nr = nr;
  rcu->mask = nr - 1lu;
  for (u64 i = 0; i < nr; i++) {
    rcu_node_init(&(rcu->nodes[i]));
  }
}

  struct rcu *
rcu_create(const u64 nr)
{
  struct rcu * const rcu = yalloc(sizeof(*rcu) + (nr * sizeof(rcu->nodes[0])));
  rcu_init(rcu, nr);
  return rcu;
}

  void *
rcu_ref(struct rcu * const rcu, const u64 magic)
{
  struct rcu_node * const node = &(rcu->nodes[magic & rcu->mask]);
  return rcu_node_ref(node);
}

  void
rcu_unref(struct rcu * const rcu, void * const ptr, const u64 magic)
{
  struct rcu_node * const node = &(rcu->nodes[magic & rcu->mask]);
  rcu_node_unref(node, ptr);
}

  void
rcu_update(struct rcu * const rcu, void * const ptr)
{
  const u64 xx = (((u64)ptr) << RCU_VALUE_SHIFT) | RCU_VALID_MASK;
  const u64 nr = rcu->nr;
  for (u64 i = 0; i < nr; i++) {
    struct rcu_node * const node = &(rcu->nodes[i]);
    for (u64 j = 0; j < 2; j++) {
      if ((node->x[j] & RCU_VALID_MASK) == 0) {
        // enable node pointer
        atomic_store(&(node->x[j]), xx);
        au64 * const v = &(node->x[1-j]);
        debug_assert(((*v) >> RCU_VALUE_SHIFT) != ((u64)ptr));
        (void)atomic_fetch_sub(v, RCU_VALID_MASK);
        break;
      }
    }
  }
  cpu_cfence();
  // wait
  for (u64 i = 0; i < nr; i++) {
    struct rcu_node * const node = &(rcu->nodes[i]);
    for (u64 j = 0; j < 2; j++) {
      const u64 x = node->x[j];
      if ((x & RCU_VALID_MASK) == 0) {
        au64 * const v = &(node->x[j]);
        while (atomic_load(v) & RCU_COUNT_MASK)
          _mm_pause();
        break;
      }
    }
  }
}
// }}} multi-rcu

// qsbr {{{
#define QSBR_STATES_NR ((492))
#define QSBR_BITMAP_NR  ((8))
struct qsbr {
  mutexlock lock;
  u64 target;
  u64 bitmap[QSBR_BITMAP_NR];
  u64 * ptrs[QSBR_STATES_NR];
};

  static inline struct qsbr *
qsbr_create(void)
{
  debug_assert(sizeof(struct qsbr) <= PGSZ);
  struct qsbr * const q = (typeof(q))pages_alloc_4kb(1);
  mutexlock_init(&(q->lock));
  return q;
}

  static inline bool
qsbr_register(struct qsbr * const q, u64 * const ptr)
{
  if (ptr == NULL) return false;
  mutexlock_lock(&(q->lock));
  for (u64 i = 0; i < QSBR_BITMAP_NR; i++) {
    const u64 bits = q->bitmap[i];
    if (bits != (~0lu)) { // bits contains zero
      const u64 pos0 = __builtin_ctzl(~bits);
      const u64 pos = (i << 6) + pos0;
      if (pos < QSBR_STATES_NR) {
        q->bitmap[i] |= (1lu << pos0);
        q->ptrs[pos] = ptr;
        mutexlock_unlock(&(q->lock));
        return true;
      }
    }
  }
  mutexlock_unlock(&(q->lock));
  return false;
}

  static inline void
qsbr_unregister(struct qsbr * const q, u64 * const ptr)
{
  while (mutexlock_trylock(&(q->lock)) == false) {
    (*ptr) = q->target;
    _mm_pause();
  }
  for (u64 i = 0; i < QSBR_BITMAP_NR; i++) {
    u64 bits = q->bitmap[i];
    while (bits) { // bits contains zero
      const u64 pos0 = __builtin_ctzl(bits);
      const u64 pos = (i << 6) + pos0;
      if (q->ptrs[pos] == ptr) {
        q->bitmap[i] &= (~(1lu << pos0));
        q->ptrs[pos] = NULL;
        mutexlock_unlock(&(q->lock));
        return;
      }
      bits &= (~(1lu << pos0));
    }
  }
  mutexlock_unlock(&(q->lock));
}

  static inline void
qsbr_wait(struct qsbr * const q, const u64 target)
{
  mutexlock_lock(&(q->lock));
  u64 bm[QSBR_BITMAP_NR] = {};
  u64 rem = 0;
  q->target = target;
  for (u64 i = 0; i < QSBR_BITMAP_NR; i++) {
    bm[i] = q->bitmap[i];
    rem += __builtin_popcountl(bm[i]);
  }

  const u64 t0 = rdtsc();
  while ((rdtsc() - t0) < 1000lu);
  // wait
  while (rem) {
    for (u64 i = 0; i < QSBR_BITMAP_NR; i++) {
      u64 bits = bm[i];
      while (bits) { // bits contains zero
        const u64 pos0 = __builtin_ctzl(bits);
        const u64 pos = (i << 6) + pos0;
        if (target == (*(q->ptrs[pos]))) {
          bm[i] &= (~(1lu << pos0));
          rem--;
        }
        bits &= (~(1lu << pos0));
      }
    }
  }
  mutexlock_unlock(&(q->lock));
}

  static inline void
qsbr_destroy(struct qsbr * const q)
{
  pages_unmap(q, PGSZ);
}
// }}} qsbr

// }}} rcu

// kv {{{
  inline size_t
kv_size(const struct kv * const kv)
{
  return sizeof(*kv) + kv->klen + kv->vlen;
}

  inline size_t
kv_size_align(const struct kv * const kv, const u64 align)
{
  debug_assert(align && ((align & (align - 1)) == 0));
  return (sizeof(*kv) + kv->klen + kv->vlen + (align - 1)) & (~(align - 1));
}

  inline size_t
key_size(const struct kv *const key)
{
  return sizeof(*key) + key->klen;
}

  inline size_t
key_size_align(const struct kv *const key, const u64 align)
{
  debug_assert(align && ((align & (align - 1)) == 0));
  return (sizeof(*key) + key->klen + (align - 1)) & (~(align - 1));
}

  inline void
kv_update_hash(struct kv * const kv)
{
  const u32 lo = crc32c((const void *)kv->kv, (const size_t)kv->klen);
  kv->hash = crc32c_extend(lo);
}

  inline void
kv_refill(struct kv * const kv, const void * const key, const u32 klen,
    const void * const value, const u32 vlen)
{
  if (kv) {
    kv->klen = klen;
    kv->vlen = vlen;
    memcpy(&(kv->kv[0]), key, klen);
    memcpy(&(kv->kv[klen]), value, vlen);
    kv_update_hash(kv);
  }
}

  inline void
kv_refill_str(struct kv * const kv, const char * const key, const char * const value)
{
  kv_refill(kv, key, (u32)strlen(key), value, (u32)strlen(value));
}

  inline struct kv *
kv_create(const void * const key, const u32 klen, const void * const value, const u32 vlen)
{
  struct kv * const kv = malloc(sizeof(*kv) + klen + vlen);
  kv_refill(kv, key, klen, value, vlen);
  return kv;
}

  inline struct kv *
kv_create_str(const char * const key, const char * const value)
{
  return kv_create(key, (u32)strlen(key), value, (u32)strlen(value));
}

  inline struct kv *
kv_dup(const struct kv * const kv)
{
  if (kv == NULL) return NULL;

  const size_t sz = kv_size(kv);
  struct kv * const new = (typeof(new))malloc(sz);
  if (new) {
    memcpy(new, kv, sz);
  }
  return new;
}

  inline struct kv *
kv_dup_key(const struct kv * const kv)
{
  if (kv == NULL) return NULL;

  const size_t sz = key_size(kv);
  struct kv * const new = (typeof(new))malloc(sz);
  if (new) {
    memcpy(new, kv, sz);
  }
  return new;
}

  inline struct kv *
kv_dup2(const struct kv * const from, struct kv * const to)
{
  if (from == NULL) return NULL;
  const size_t sz = kv_size(from);
  struct kv * const new = to ? to : (typeof(new))malloc(sz);
  memcpy(new, from, sz);
  return new;
}

  inline struct kv *
kv_dup2_key(const struct kv * const from, struct kv * const to)
{
  if (from == NULL) return NULL;
  const size_t sz = key_size(from);
  struct kv * const new = to ? to : (typeof(new))malloc(sz);
  memcpy(new, from, sz);
  new->vlen = 0;
  return new;
}

  inline struct kv *
kv_dup2_key_prefix(const struct kv * const from, struct kv * const to, const u64 plen)
{
  if (from == NULL) return NULL;
  const size_t sz = key_size(from) - from->klen + plen;
  struct kv * const new = to ? to : (typeof(new))malloc(sz);
  if (new) {
    new->klen = plen;
    memcpy(new->kv, from->kv, plen);
    new->vlen = 0;
    kv_update_hash(new);
  }
  return new;
}

  inline struct kv *
kv_alloc_malloc(const u64 size, void * const priv)
{
  (void)priv;
  return (struct kv *)malloc(size);
}

  void
kv_retire_free(struct kv * const kv, void * const priv)
{
  (void)priv;
  free(kv);
}

// key1 and key2 must be valid ptr
  inline bool
kv_keymatch(const struct kv * const key1, const struct kv * const key2)
{
  if ((key1->hash == key2->hash) && (key1->klen == key2->klen) && (0 == memcmp(key1->kv, key2->kv, key1->klen)))
    return true;
  return false;
}

  inline bool
kv_fullmatch(const struct kv * const kv1, const struct kv * const kv2)
{
  if ((kv1->kvlen == kv2->kvlen) && (0 == memcmp(kv1, kv2, sizeof(*kv1) + kv1->klen + kv1->vlen)))
    return true;
  return false;
}

  int
kv_keycompare(const struct kv * const kv1, const struct kv * const kv2)
{
  debug_assert(kv1);
  debug_assert(kv2);
  const u32 len = kv1->klen < kv2->klen ? kv1->klen : kv2->klen;
  const int cmp = memcmp(kv1->kv, kv2->kv, (size_t)len);
  if (cmp == 0) {
    if (kv1->klen < kv2->klen) {
      return -1;
    } else if (kv1->klen > kv2->klen) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return cmp;
  }
}

  static int
__kv_compare_pp(const void * const p1, const void * const p2)
{
  const struct kv ** const pp1 = (typeof(pp1))p1;
  const struct kv ** const pp2 = (typeof(pp2))p2;
  return kv_keycompare(*pp1, *pp2);
}

  inline void
kv_qsort(const struct kv ** const kvs, const size_t nr)
{
  qsort(kvs, nr, sizeof(kvs[0]), __kv_compare_pp);
}

  inline void *
kv_value_ptr(struct kv * const kv)
{
  return (void *)(&(kv->kv[kv->klen]));
}

  inline void *
kv_key_ptr(struct kv * const kv)
{
  return (void *)(&(kv->kv[0]));
}

  inline const void *
kv_value_ptr_const(const struct kv * const kv)
{
  return (const void *)(&(kv->kv[kv->klen]));
}

  inline const void *
kv_key_ptr_const(const struct kv * const kv)
{
  return (const void *)(&(kv->kv[0]));
}

  static inline u32
__lcp16_aligned(const void * const p1, const void * const p2)
{
  const __m128i v1 = _mm_load_si128(p1);
  const __m128i v2 = _mm_load_si128(p2);
  const u32 mask = ~((u32)_mm_movemask_epi8(_mm_cmpeq_epi8(v1, v2)));
  return __builtin_ctz(mask);
}

// return the length of longest common prefix of the two keys
  inline u32
kv_key_lcp(const struct kv * const key1, const struct kv * const key2)
{
  const u32 max = (key1->klen < key2->klen) ? key1->klen : key2->klen;
  const u32 max128 = max & (~0xfu);
  u32 clen = 0;
  const u8 * p1 = key1->kv;
  const u8 * p2 = key2->kv;
  // inc by 4
  while (clen < max128) {
    const u32 lcpinc = __lcp16_aligned(p1, p2);
    if (lcpinc < 16) return clen + lcpinc;
    clen += 16;
    p1 += 16;
    p2 += 16;
  }
  const u32 max32 = max & (~0x3u);
  // inc by 4
  while (clen < max32) {
    const u32 v1 = *(const u32 *)p1;
    const u32 v2 = *(const u32 *)p2;
    if (v1 != v2) return clen + (__builtin_ctz(v1 ^ v2) >> 3);
    clen += sizeof(u32);
    p1 += sizeof(u32);
    p2 += sizeof(u32);
  }
  while ((clen < max) && (*p1 == *p2)) {
    clen++;
    p1++;
    p2++;
  }
  return clen;
}

// return true if p is a prefix of key
  inline bool
kv_key_is_prefix(const struct kv * const p, const struct kv * const key)
{
  return ((p->klen <= key->klen) && (kv_key_lcp(p, key) == p->klen)) ? true : false;
}

  static inline struct kv *
u64_to_kv(const u64 v)
{
  return (struct kv *)v;
}

  static inline u64
kv_to_u64(const struct kv * const ptr)
{
  return (u64)ptr;
}

  static inline u64
kvmap_pkey(const u64 hash)
{
  return ((hash >> 16) ^ hash) & UINT64_C(0xffff);
}

  static const char *
__kv_pattern(const char c)
{
  switch (c) {
    case 's': return "%c";
    case 'x': return " %02hhx";
    case 'd': return " %03hhu";
    default: return NULL;
  }
}

// cmd "KV" K and V can be 's' for string, 'x' for hex, 'd' for dec.
// n for newline after kv
  void
kv_print(const struct kv * const kv, const char * const cmd, FILE * const out)
{
  debug_assert(cmd);
  const u32 klen = kv->klen;
  fprintf(out, "#%04lx #%016lx k[%2u] ", kvmap_pkey(kv->hash), kv->hash, klen);
  const u32 klim = klen < 1024u ? klen : 1024u;

  const char * const kpat = __kv_pattern(cmd[0]);
  for (u32 i = 0; i < klim; i++)
    fprintf(out, kpat, kv->kv[i]);
  if (klim < klen)
    fprintf(out, " ...");

  const char * const vpat = __kv_pattern(cmd[1]);
  if (vpat) {
    const u32 vlen = kv->vlen;
    const u32 vlim = vlen < 1024u ? vlen : 1024u;
    fprintf(out, "  v[%4u] ", vlen);
    for (u32 i = 0; i < vlim; i++)
      fprintf(out, vpat, kv->kv[klen + i]);
    if (vlim < vlen)
      fprintf(out, " ...");
  }
  if (strchr(cmd, 'n'))
    fprintf(out, "\n");
}
// }}} kv

// kvmap {{{
struct kv_entry {
  union {
    struct {
      u64 pkey:16;
      u64 ptr:48;
    };
    u64 v64;
  };
};
#define KVMAP_SLOT_NR ((8lu))
struct kvmap_slot {
  struct kv_entry e[KVMAP_SLOT_NR];
};

#define KVMAP_BASE_LEVEL ((6))
#define KVMAP_HLOCK_POWER ((10))
#define KVMAP_HLOCK_NR ((UINT64_C(1) << KVMAP_HLOCK_POWER))
#define KVMAP_HLOCK_MASK ((KVMAP_HLOCK_NR - 1))
#define KVMAP_MIN_LEVEL ((KVMAP_BASE_LEVEL + KVMAP_HLOCK_POWER))
#define KVMAP_SLOTS_OF(__lid__) ((UINT64_C(1) << (__lid__)))
#define KVMAP_ENTRIES_OF(__lid__) ((UINT64_C(8) << (__lid__)))
#define KVMAP_PAGES_OF(__lid__) ((UINT64_C(1) << ((__lid__) - UINT64_C(6))))
#define KVMAP_SIZE_OF_LEVEL(__lid__) ((KVMAP_PAGES_OF(__lid__) * PGSZ))

  static int
__kvmap_entry_pkey_comp(const void * const p1, const void * const p2)
{
  const struct kv_entry * const e1 = (typeof(e1))p1;
  const struct kv_entry * const e2 = (typeof(e2))p2;
  if (e1->pkey < e2->pkey) {
    return -1;
  } else if (e1->pkey > e2->pkey) {
    return 1;
  } else {
    return 0;
  }
}

  static inline int
__kvmap_entry_keycompare_vptr(const void * const p1, const void * const p2)
{
  const struct kv_entry * const e1 = (typeof(e1))p1;
  const struct kv_entry * const e2 = (typeof(e2))p2;
  const struct kv * const k1 = u64_to_kv(e1->ptr);
  const struct kv * const k2 = u64_to_kv(e2->ptr);
  return kv_keycompare(k1, k2);
}

  static inline void
__kvmap_entry_qsort(struct kv_entry * const es, const size_t nr)
{
  qsort(es, nr, sizeof(es[0]), __kvmap_entry_keycompare_vptr);
}

  static inline void
kvmap_put_entry(struct kvmap_mm * const mm, struct kv_entry * const e, const struct kv * const kv)
{
  struct kv * const old = u64_to_kv(e->ptr);
  if (old && mm->rf)
    mm->rf(old, mm->rp);
  if (kv) {
    e->ptr = kv_to_u64(kv);
    e->pkey = kvmap_pkey(kv->hash);
  } else {
    e->v64 = 0lu;
  }
}

static const struct kvmap_mm __kvmap_mm_default = {
  .af = kv_alloc_malloc,
  .ap = NULL,
  .rf = kv_retire_free,
  .rp = NULL,
};
// }}} kvmap

// wormhole {{{

// def {{{
#define WORMHOLE_BASE ((16u))
#define WORMHOLE_KPN ((WORMHOLE_BASE * 8u)) // keys per node; power of 2
#define WORMHOLE_MID ((WORMHOLE_BASE * 4u)) // ideal cut point for split, the closer the better

#define WORMHOLE_KPN_MIN ((WORMHOLE_BASE * 2)) // 2/8
#define WORMHOLE_KPN_MRG ((WORMHOLE_BASE * 6)) // 6/8

// FO is fixed at 256. Don't change it
#define WORMHOLE_FO  ((256)) // index fan-out
// number of bits in a bitmap
#define WORMHOLE_BMNR ((WORMHOLE_FO >> 6)) // number of u64
// }}} def

// struct {{{
struct wormmeta {
  struct kv * key;
  u32 bitcount;
  u32 bitmin;
  struct wormleaf * lmost;
  struct wormleaf * rmost;
  u64 bitmap[WORMHOLE_BMNR];
};

struct wormleaf {
  // first line
  struct wormleaf * prev; // prev leaf
  struct wormleaf * next; // next leaf
  struct kv * anchor;
  u32 nr_sorted;
  u32 nr_keys;
  u64 version;
  u32 klen;
  u32 padding0;
  u64 padding1;
  // second line
  rwlock leaflock;
  struct kv_entry eh[WORMHOLE_KPN]; // sorted by hashes
  struct kv_entry es[WORMHOLE_KPN]; // sorted by keys
};

struct wormslot {
  u16 t[KVMAP_SLOT_NR];
}__attribute__((packed));

struct wormhmap {
  u64 version;
  struct wormslot * wmap;
  u64 mask;
  struct kvmap_slot * pmap;

  u32 maxplen;
  u32 stale;
  struct wormhole * map;
  u64 psize;
  u64 wsize;
};

struct wormhole {
  // 1 line
  struct wormhmap * hmap;
  u64 padding0[7];
  // 1 line
  struct kvmap_mm mm;
  struct wormleaf * leaf0;
  struct qsbr * qsbr;
  struct slab * slab_meta;
  struct slab * slab_leaf;
  // 2 lines
  struct wormhmap hmap2[2];
  // fifth line
  rwlock metalock;
  u64 padding1[7];
};

struct wormhole_iter {
  struct wormref * ref;
  struct wormleaf * leaf;
  u32 next_id;
};

struct wormkref { // reference to a key
  u64 hash;
  u32 plen; // prefix length; plen <= klen
  u32 klen; // the original klen
  const u8 * key;
};

struct wormref {
  struct wormhole * map;
  u64 qstate;
};
// }}} struct

// helpers {{{

// alloc {{{
  static inline struct kv *
wormhole_alloc_key(const size_t klen)
{
  // evaluation says slab and yalloc are worse...
  return (struct kv *)malloc(sizeof(struct kv) + klen);
}

  static inline struct kv *
wormhole_alloc_kv(struct wormhole * const map, const size_t klen, const size_t vlen)
{
  const size_t size = sizeof(struct kv) + klen + vlen;
  return map->mm.af(size, map->mm.ap);
}

  static inline void
wormhole_free_key(struct kv * const key)
{
  free(key);
}

  static struct wormleaf *
wormhole_alloc_leaf(struct wormhole * const map, struct wormleaf * const prev,
    struct wormleaf * const next, struct kv * const anchor)
{
  struct wormleaf * const leaf = (typeof(leaf))slab_alloc(map->slab_leaf);
  debug_assert(leaf);
  rwlock_init(&(leaf->leaflock));
  leaf->version = 0;
  leaf->anchor = anchor;
  if (anchor)
    leaf->klen = anchor->klen;
  leaf->nr_sorted = 0u;
  leaf->nr_keys = 0u;
  leaf->prev = prev;
  leaf->next = next;
  return leaf;
}

  static inline struct wormmeta *
wormhole_alloc_meta(struct wormhole * const map, struct wormleaf * const lrmost,
    struct kv * const key)
{
  struct wormmeta * const meta = (typeof(meta))slab_alloc(map->slab_meta);
  debug_assert(meta);
  meta->key = key;
  meta->bitcount = 0;
  meta->bitmin = WORMHOLE_FO;
  meta->lmost = lrmost;
  meta->rmost = lrmost;
  for (u64 i = 0; i < WORMHOLE_BMNR; i++)
    meta->bitmap[i] = 0;
  return meta;
}
// }}} alloc

// key/prefix {{{
  static inline bool
wormhole_kref_keymatch(const struct wormkref * const ref, const struct kv * const key)
{
  if ((ref->hash == key->hash) && (ref->plen == key->klen)) {
    return memcmp(ref->key, key->kv, ref->plen) ? false : true;
  }
  return false;
}

// warning: buffer overflow risk, call it carefully
  static inline void
wormhole_prefix(struct kv * const prefix, const u32 klen)
{
  prefix->klen = klen;
  kv_update_hash(prefix);
}

  static inline void
wormhole_prefix_inc_short(struct kv * const prefix, const u32 klen)
{
  debug_assert(klen >= prefix->klen);
  const u32 lo = crc32c_inc_short(prefix->kv + prefix->klen, klen - prefix->klen, (u32)prefix->hash);
  prefix->hash = crc32c_extend(lo);
  prefix->klen = klen;
}

// this is for meta_up
  static inline void
wormhole_kref_inc_long(struct wormkref * const ref, const u32 plen, const u32 seed, const u32 slen)
{
  ref->hash = crc32c_extend(crc32c_inc(ref->key + slen, plen - slen, seed));
  ref->plen = plen;
}

  static inline void
wormhole_kref_inc_short(struct wormkref * const ref, const u32 plen, const u32 seed, const u32 slen)
{
  ref->hash = crc32c_extend(crc32c_inc_short(ref->key + slen, plen - slen, seed));
  ref->plen = plen;
}
// }}} key/prefix

// meta/bitmap {{{
  static inline bool
wormhole_meta_bm_test(const struct wormmeta * const meta, const u32 id)
{
  return (bool)((meta->bitmap[id >> 6] >> (id & 0x3fu)) & 1lu);
}

  static inline void
wormhole_meta_bm_set(struct wormmeta * const meta, const u32 id)
{
  meta->bitmap[id >> 6u] |= (1lu << (id & 0x3fu));
  meta->bitcount = _mm_popcnt_u64(meta->bitmap[0]) + _mm_popcnt_u64(meta->bitmap[1])
    + _mm_popcnt_u64(meta->bitmap[2]) + _mm_popcnt_u64(meta->bitmap[3]);
  if (id < meta->bitmin)
    meta->bitmin = id;
}

  static inline u32
wormhole_meta_bm_gt(const struct wormmeta * const meta, const u32 id0)
{
  if ((id0 & 0x3fu) != 0x3fu) { // not at 63
    const u32 id = id0 + 1u;
    const u64 bits = meta->bitmap[id >> 6] >> (id & 0x3fu);
    if (bits)
      return id + (u32)__builtin_ctzl(bits);
  }
  for (u32 ix = (id0 >> 6) + 1; ix < 4; ix++) {
    if (meta->bitmap[ix])
      return (ix << 6) + (u32)(__builtin_ctzl(meta->bitmap[ix]));
  }
  return WORMHOLE_FO;
}

  static inline void
wormhole_meta_bm_clear(struct wormmeta * const meta, const u32 id)
{
  meta->bitmap[id >> 6u] &= (~(1lu << (id & 0x3fu)));
  meta->bitcount = _mm_popcnt_u64(meta->bitmap[0]) + _mm_popcnt_u64(meta->bitmap[1])
    + _mm_popcnt_u64(meta->bitmap[2]) + _mm_popcnt_u64(meta->bitmap[3]);
  if (id == meta->bitmin) {
    meta->bitmin = wormhole_meta_bm_gt(meta, id);
    debug_assert(meta->bitmin > id);
  }
}

  static inline struct wormmeta *
__u64_to_meta(const u64 v)
{
  return (struct wormmeta *)v;
}

  static inline u64
__meta_to_u64(const struct wormmeta * const ptr)
{
  return (u64)ptr;
}

// find the highest bit that is lower to the given id in the bitmap
// returns original id if not found
  static inline u32
wormhole_meta_bm_lt(const struct wormmeta * const meta, const u32 id0)
{
  if (id0 & 0x3fu) { // not at 0
    const u32 id = id0 - 1u;
    const u64 bits = meta->bitmap[id >> 6] << (63u - (id & 0x3fu));
    if (bits)
      return id - (u32)__builtin_clzl(bits);
  }
  for (u32 ixp = id0 >> 6; ixp; ixp--) {
    if (meta->bitmap[ixp-1u])
      return (ixp << 6) - 1u - (u32)(__builtin_clzl(meta->bitmap[ixp-1u]));
  }
  return id0;
}
// }}} meta/bitmap

// }}} helpers

// hmap {{{
  static inline u16
wormhole_hmap_skey(const u16 pkey)
{
  return pkey | (!pkey);
  //return pkey | 1;
}

  static inline u64
wormhole_hmap_peek_slot(const struct wormslot * const s, const __m128i skey)
{
  return (u64)_mm_movemask_epi8(_mm_cmpeq_epi16(skey, _mm_load_si128((const void *)(s))));
}

  static inline bool
wormhole_hmap_peek(const struct wormhmap * const hmap, const u64 hash)
{
  const __m128i sk = _mm_set1_epi16(wormhole_hmap_skey(kvmap_pkey(hash)));
  const u64 midx = hash & hmap->mask;
  const u64 midy = bswap_64(hash) & hmap->mask;
  if (wormhole_hmap_peek_slot(&(hmap->wmap[midx]), sk))
    return true;
  return wormhole_hmap_peek_slot(&(hmap->wmap[midy]), sk);
}

  static inline u64
wormhole_hmap_count_entry(const struct wormhmap * const hmap, const u64 mid)
{
  const __m128i skey = _mm_set1_epi16((u16)0u);
  const u64 mask = wormhole_hmap_peek_slot(&(hmap->wmap[mid]), skey);
  return mask ? (__builtin_ctzl(mask) >> 1) : 8;
}

  static inline struct wormmeta *
wormhole_hmap_get_slot(const struct wormhmap * const hmap, const u64 mid, const __m128i skey,
    const struct kv * const key)
{
  const struct kvmap_slot * const pmap = hmap->pmap;
  u64 mask = wormhole_hmap_peek_slot(&(hmap->wmap[mid]), skey);
  while (mask) {
    const u32 i = __builtin_ctzl(mask) >> 1;
    struct wormmeta * const meta = __u64_to_meta(pmap[mid].e[i].ptr);
    if (kv_keymatch(key, meta->key))
      return meta;
    mask ^= (3lu << (i << 1));
  }
  return NULL;
}

  static inline struct wormmeta *
wormhole_hmap_get(const struct wormhmap * const hmap, const struct kv * const key)
{
  const u64 hash = key->hash;
  const u64 midx = hash & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midx]), 0, 0);
  const u64 midy = bswap_64(hash) & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midy]), 0, 0);
  const __m128i skey = _mm_set1_epi16(wormhole_hmap_skey(kvmap_pkey(hash)));

  struct wormmeta * const r = wormhole_hmap_get_slot(hmap, midx, skey, key);
  if (r)
    return r;
  return wormhole_hmap_get_slot(hmap, midy, skey, key);
}

  static inline struct wormmeta *
wormhole_hmap_get_kref_slot(const struct wormhmap * const hmap, const u64 mid, const __m128i skey,
    const struct wormkref * const ref)
{
  const struct kvmap_slot * const pmap = hmap->pmap;
  u64 mask = wormhole_hmap_peek_slot(&(hmap->wmap[mid]), skey);
  while (mask) {
    const u32 i = __builtin_ctzl(mask) >> 1;
    struct wormmeta * const meta = __u64_to_meta(pmap[mid].e[i].ptr);
    if (wormhole_kref_keymatch(ref, meta->key)) {
      if ((ref->klen == ref->plen) || (meta->bitmin > ref->key[ref->plen])) {
        __builtin_prefetch(meta->lmost, 0, 0);
      }
      return meta;
    }
    mask ^= (3lu << (i << 1));
  }
  return NULL;
}

  static inline struct wormmeta *
wormhole_hmap_get_kref(const struct wormhmap * const hmap, const struct wormkref * const ref)
{
  const u64 hash = ref->hash;
  const u64 midx = hash & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midx]), 0, 0);
  const u64 midy = bswap_64(hash) & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midy]), 0, 0);
  const __m128i skey = _mm_set1_epi16(wormhole_hmap_skey(kvmap_pkey(hash)));

  struct wormmeta * const r = wormhole_hmap_get_kref_slot(hmap, midx, skey, ref);
  if (r)
    return r;
  return wormhole_hmap_get_kref_slot(hmap, midy, skey, ref);
}

  static inline struct wormmeta *
wormhole_hmap_get_kref1_slot(const struct wormhmap * const hmap, const u64 mid, const __m128i skey,
    const struct wormkref * const ref, const u32 cid)
{
  const struct kvmap_slot * const pmap = hmap->pmap;
  const u32 plen = ref->plen;
  u64 mask = wormhole_hmap_peek_slot(&(hmap->wmap[mid]), skey);
  while (mask) {
    const u32 i = __builtin_ctzl(mask) >> 1;
    struct wormmeta * const meta = __u64_to_meta(pmap[mid].e[i].ptr);
    const struct kv * const mkey = meta->key;
    if (mkey->klen == (plen + 1u) && (0 == memcmp(ref->key, mkey->kv, plen))
        && (mkey->kv[plen] == cid)) {
      __builtin_prefetch(meta->rmost, 0, 0);
      return meta;
    }
    mask ^= (3lu << (i << 1));
  }
  return NULL;
}

  static inline struct wormmeta *
wormhole_hmap_get_kref1(const struct wormhmap * const hmap, const struct wormkref * const ref, const u32 cid)
{
  const u64 hash = crc32c_extend(_mm_crc32_u8((u32)ref->hash, cid));
  const u64 midx = hash & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midx]), 0, 0);
  const u64 midy = bswap_64(hash) & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midy]), 0, 0);
  const __m128i skey = _mm_set1_epi16(wormhole_hmap_skey(kvmap_pkey(hash)));

  struct wormmeta * const r = wormhole_hmap_get_kref1_slot(hmap, midx, skey, ref, cid);
  if (r)
    return r;
  return wormhole_hmap_get_kref1_slot(hmap, midy, skey, ref, cid);
}

  static inline void
wormhole_hmap_insert_at(const struct wormhmap * const hmap, const u64 mid,
    const struct wormmeta * const meta, const u64 pkey, const u64 i)
{
  debug_assert(hmap->pmap[mid].e[i].v64 == 0lu);
  debug_assert(hmap->wmap[mid].t[i] == 0u);
  hmap->pmap[mid].e[i].pkey = pkey;
  hmap->pmap[mid].e[i].ptr = __meta_to_u64(meta);
  hmap->wmap[mid].t[i] = wormhole_hmap_skey(pkey);
}

  static inline void
wormhole_hmap_insert(const struct wormhmap * const hmap, const u64 mid,
    const struct wormmeta * const meta, const u64 pkey)
{
  const u64 i = wormhole_hmap_count_entry(hmap, mid);
  debug_assert(i < 8lu);
  wormhole_hmap_insert_at(hmap, mid, meta, pkey, i);
}

  static void
wormhole_hmap_squeeze(const struct wormhmap * const hmap)
{
  const u64 nrs = hmap->mask + 1lu;
  struct kvmap_slot * const pmap = hmap->pmap;
  struct wormslot * const wmap = hmap->wmap;
  const u64 mask = hmap->mask;
  for (u64 si = 0; si < nrs; si++) {
    for (u64 i = 0; i < KVMAP_SLOT_NR; i++) {
      const u64 ei = KVMAP_SLOT_NR - i - 1lu;
      if (pmap[si].e[ei].v64) {
        struct wormmeta * const meta = __u64_to_meta(pmap[si].e[ei].ptr);
        const u64 sj = meta->key->hash & mask;
        if (sj != si) {
          const u64 ej = wormhole_hmap_count_entry(hmap, sj);
          if (ej < KVMAP_SLOT_NR) {
            pmap[sj].e[ej] = pmap[si].e[ei];
            wmap[sj].t[ej] = wmap[si].t[ei];
            const u64 ni = wormhole_hmap_count_entry(hmap, si) - 1lu;
            pmap[si].e[ei] = pmap[si].e[ni];
            wmap[si].t[ei] = wmap[si].t[ni];
            pmap[si].e[ni].v64 = 0lu;
            wmap[si].t[ni] = 0u;
          }
        }
      }
    }
  }
}

  static bool
wormhole_hmap_expand(struct wormhmap * const hmap)
{
  // sync expand
  const u64 mask0 = hmap->mask;
  const u64 nr0 = mask0 + 1lu;
  const u64 mask1 = mask0 + nr0;
  const u64 nr1 = nr0 << 1;
  u64 psize = nr1 * sizeof(hmap->pmap[0]);
  u64 wsize = nr1 * sizeof(hmap->wmap[0]);
  struct kvmap_slot * const ps = pages_alloc_best(psize, true, &psize);
  struct wormslot * const ws = pages_alloc_best(wsize, true, &wsize);
  if ((ps == NULL) || (ws == NULL)) {
    if (ps)
      pages_unmap(ps, psize);
    if (ws)
      pages_unmap(ws, wsize);
    return false;
  }

  struct wormhmap hmap1 = *hmap;
  hmap1.pmap = ps;
  hmap1.psize = psize;
  hmap1.wmap = ws;
  hmap1.wsize = wsize;
  hmap1.mask = mask1;

  const struct kvmap_slot * const pmap0 = hmap->pmap;

  for (u64 s = 0; s < nr0; s++) {
    const struct kv_entry * e = &(pmap0[s].e[0]);
    for (u64 i = 0; (i < KVMAP_SLOT_NR) && e->v64; i++, e++) {
      const struct wormmeta * const meta = __u64_to_meta(e->ptr);
      const struct kv * const mkey = meta->key;
      const u64 pkey = kvmap_pkey(mkey->hash);
      const u64 idx0 = mkey->hash & mask0;
      const u64 idx1 = ((idx0 == s) ? mkey->hash : bswap_64(mkey->hash)) & mask1;
      wormhole_hmap_insert(&hmap1, idx1, meta, pkey);
    }
  }
  pages_unmap(hmap->pmap, hmap->psize);
  pages_unmap(hmap->wmap, hmap->wsize);
  *hmap = hmap1;
  wormhole_hmap_squeeze(hmap);
  return true;
}

  static bool
wormhole_hmap_update_slot(struct wormhmap * const hmap, const u64 mid, const struct kv * const key,
    const struct wormmeta * const meta, const __m128i skey)
{
  u64 mask = wormhole_hmap_peek_slot(&(hmap->wmap[mid]), skey);
  while (mask) {
    const u32 i = __builtin_ctzl(mask) >> 1;
    const struct wormmeta * const meta0 = __u64_to_meta(hmap->pmap[mid].e[i].ptr);
    if (kv_keymatch(key, meta0->key)) {
      hmap->pmap[mid].e[i].ptr = __meta_to_u64(meta);
      return true;
    }
    mask -= (3u << (i << 1));
  }
  return false;
}

  static inline u64
wormhole_hmap_cuckoo(struct wormhmap * const hmap, const u64 mid0, const u64 depth)
{
  struct kv_entry * e = &(hmap->pmap[mid0].e[0]);
  for (u64 i = 0; (i < KVMAP_SLOT_NR) && e->v64; i++, e++) {
    const struct wormmeta * const meta = __u64_to_meta(e->ptr);
    const struct kv * const mkey = meta->key;

    const u64 midx = mkey->hash & hmap->mask;
    const u64 midy = bswap_64(mkey->hash) & hmap->mask;
    const u64 midt = (midx != mid0) ? midx : midy;
    if (midt != mid0) {
      const u64 c1 = wormhole_hmap_count_entry(hmap, midt);
      const u64 ct = ((c1 == KVMAP_SLOT_NR) && depth) ? wormhole_hmap_cuckoo(hmap, midt, depth - 1) : c1;
      if (ct < KVMAP_SLOT_NR) {
        hmap->pmap[midt].e[ct] = *e;
        hmap->wmap[midt].t[ct] = wormhole_hmap_skey(e->pkey);
        return i;
      }
    }
  }
  return KVMAP_SLOT_NR;
}

  static bool
wormhole_hmap_set(struct wormhmap * const hmap, const struct kv * const key,
    const struct wormmeta * const meta)
{
  const u64 hash = key->hash;
  const u64 midx = hash & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midx]), 0, 0);
  const u64 midy = bswap_64(hash) & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midy]), 0, 0);
  const u64 pkey = kvmap_pkey(hash);
  const __m128i skey = _mm_set1_epi16(wormhole_hmap_skey(pkey));
  if (wormhole_hmap_update_slot(hmap, midx, key, meta, skey))
    return true;

  if (wormhole_hmap_update_slot(hmap, midy, key, meta, skey))
    return true;

  // insert
  const u64 i1 = wormhole_hmap_count_entry(hmap, midx);
  if (i1 < KVMAP_SLOT_NR) {
    wormhole_hmap_insert_at(hmap, midx, meta, pkey, i1);
    return true;
  }

  const u64 i2 = wormhole_hmap_count_entry(hmap, midy);
  if (i2 < KVMAP_SLOT_NR) {
    wormhole_hmap_insert_at(hmap, midy, meta, pkey, i2);
    return true;
  }

  // cuckoo
  const u64 c1 = wormhole_hmap_cuckoo(hmap, midx, 1);
  if (c1 < KVMAP_SLOT_NR) {
    struct kvmap_slot * const s1 = &(hmap->pmap[midx]);
    s1->e[c1].ptr = __meta_to_u64(meta);
    s1->e[c1].pkey = pkey;
    hmap->wmap[midx].t[c1] = wormhole_hmap_skey(pkey);
    return true;
  }
  const u64 c2 = wormhole_hmap_cuckoo(hmap, midy, 1);
  if (c2 < KVMAP_SLOT_NR) {
    struct kvmap_slot * const s2 = &(hmap->pmap[midy]);
    s2->e[c2].ptr = __meta_to_u64(meta);
    s2->e[c2].pkey = pkey;
    hmap->wmap[midy].t[c2] = wormhole_hmap_skey(pkey);
    return true;
  }
  // expand
  if (wormhole_hmap_expand(hmap) == false)
    return false;

  return wormhole_hmap_set(hmap, key, meta);
}

  static bool
wormhole_hmap_del_slot(struct wormhmap * const hmap, const u64 mid,
    const struct kv * const key, const __m128i skey)
{
  u64 mask = wormhole_hmap_peek_slot(&(hmap->wmap[mid]), skey);
  while (mask) {
    const u32 i = __builtin_ctzl(mask) >> 1;
    const struct wormmeta * const meta = __u64_to_meta(hmap->pmap[mid].e[i].ptr);
    if (kv_keymatch(key, meta->key)) {
      const u64 j = wormhole_hmap_count_entry(hmap, mid) - 1lu;
      hmap->pmap[mid].e[i] = hmap->pmap[mid].e[j];
      hmap->pmap[mid].e[j].v64 = 0lu;
      hmap->wmap[mid].t[i] = hmap->wmap[mid].t[j];
      hmap->wmap[mid].t[j] = 0u;
      return true;
    }
    mask -= (3u << (i << 1));
  }
  return false;
}

  static bool
wormhole_hmap_del(struct wormhmap * const hmap, const struct kv * const key)
{
  const u64 hash = key->hash;
  const u64 midx = hash & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midx]), 0, 0);
  const u64 midy = bswap_64(hash) & hmap->mask;
  __builtin_prefetch(&(hmap->pmap[midy]), 0, 0);
  const __m128i skey = _mm_set1_epi16(wormhole_hmap_skey(kvmap_pkey(hash)));
  return wormhole_hmap_del_slot(hmap, midx, key, skey)
    || wormhole_hmap_del_slot(hmap, midy, key, skey);
}

  static void
wormhole_hmap_init(struct wormhmap * const hmap, struct wormhole * const map)
{
  const u64 nr = 1lu << KVMAP_MIN_LEVEL;
  const u64 psize = sizeof(hmap->pmap[0]) * nr;
  const u64 wsize = sizeof(hmap->wmap[0]) * nr;
  hmap->pmap = pages_alloc_best(psize, true, &(hmap->psize));
  hmap->wmap = pages_alloc_best(wsize, true, &(hmap->wsize));
  debug_assert(hmap->pmap);
  debug_assert(hmap->wmap);
  hmap->mask = nr - 1lu;
  hmap->version = 0;
  hmap->map = map;
  hmap->maxplen = 0u;
  hmap->stale = 0u;
}
// }}} hmap

// create {{{
  static void
wormhole_create_leaf0(struct wormhole * const map)
{
  // create leaf of empty key
  struct kv * const anchor = wormhole_alloc_key(0);
  kv_refill(anchor, NULL, 0, NULL, 0);
  struct wormleaf * const leaf0 = wormhole_alloc_leaf(map, NULL, NULL, anchor);
  map->leaf0 = leaf0;

  // create meta of empty key
  for (u64 i = 0; i < 2; i++) {
    struct kv * const mkey = wormhole_alloc_key(0);
    kv_refill(mkey, NULL, 0, NULL, 0);
    struct wormmeta * const meta0 = wormhole_alloc_meta(map, leaf0, mkey);
    const bool rset = wormhole_hmap_set(&(map->hmap2[i]), mkey, meta0);
    (void)rset;
    debug_assert(rset);
  }
}

  struct wormhole *
wormhole_create(const struct kvmap_mm * const mm)
{
  struct wormhole * const map = (typeof(map))yalloc(sizeof(*map));
  memset(map, 0, sizeof(*map));
  // mm
  map->mm = mm ? (*mm) : __kvmap_mm_default;

  // hmap
  for (u64 i = 0; i < 2; i++)
    wormhole_hmap_init(&(map->hmap2[i]), map);

  // slabs
  map->slab_meta = slab_create(sizeof(struct wormmeta), 1lu << 21);
  map->slab_leaf = slab_create(sizeof(struct wormleaf), 1lu << 24);

  // others
  map->qsbr = qsbr_create();
  wormhole_create_leaf0(map);
  rwlock_init(&(map->metalock));
  map->hmap = &(map->hmap2[0]);
  return map;
}
// }}} create

// get {{{

// jump {{{
// search for a wormmeta in the hash table that has the longest prefix of the requested key
// the returned wormmeta could be either a meta or a leaf
// the corresponding prefix is left at [pbuf] for callers to use at return
  static struct wormmeta *
wormhole_meta_up(const struct wormhmap * const hmap, struct wormkref * const ref)
{
  // invariant: lo <= lp < hi
  // finish condition: (lo + 1) == hi
  u64 lo = 0;
  u64 hi = (hmap->maxplen < ref->klen ? hmap->maxplen : ref->klen) + 1;
  u64 seed = crc32c_extend(CRC32C_SEED);

#define META_UP_GAP_1 ((7))
  while ((lo + META_UP_GAP_1) < hi) {
    const u64 pm = (lo + hi) >> 1;
    wormhole_kref_inc_long(ref, pm, (u32)seed, lo);
    if (wormhole_hmap_peek(hmap, ref->hash)) {
      seed = ref->hash;
      lo = pm;
    } else {
      hi = pm;
    }
  }

  while ((lo + 1) < hi) {
    const u64 pm = (lo + hi) >> 1;
    wormhole_kref_inc_short(ref, pm, (u32)seed, lo);
    if (wormhole_hmap_peek(hmap, ref->hash)) {
      seed = ref->hash;
      lo = pm;
    } else {
      hi = pm;
    }
  }

  if (ref->plen != lo) {
    ref->hash = seed;
    ref->plen = lo;
  }
  struct wormmeta * ret = wormhole_hmap_get_kref(hmap, ref);
  if (ret)
    return ret;

  hi = lo;
  lo = 0lu;
  seed = crc32c_extend(CRC32C_SEED);

#define META_UP_GAP_2 ((5))
  while ((lo + META_UP_GAP_2) < hi) {
    const u64 pm = (lo + hi + hi + hi) >> 2;
    wormhole_kref_inc_long(ref, pm, (u32)seed, lo);
    struct wormmeta * const tmp = wormhole_hmap_get_kref(hmap, ref);
    if (tmp) {
      seed = ref->hash;
      lo = pm;
      ret = tmp;
      if (wormhole_meta_bm_test(tmp, ref->key[pm])) {
        lo++;
        seed = crc32c_extend(_mm_crc32_u8((u32)seed, ref->key[pm]));
        ret = NULL;
      } else {
        hi = pm + 1;
        break;
      }
    } else {
      hi = pm;
    }
  }

  while ((lo + 1) < hi) {
    const u64 pm = (lo + hi + hi + hi) >> 2;
    wormhole_kref_inc_short(ref, pm, (u32)seed, lo);
    struct wormmeta * const tmp = wormhole_hmap_get_kref(hmap, ref);
    if (tmp) {
      seed = ref->hash;
      lo = pm;
      ret = tmp;
      if (wormhole_meta_bm_test(tmp, ref->key[pm])) {
        lo++;
        seed = crc32c_extend(_mm_crc32_u8((u32)seed, ref->key[pm]));
        ret = NULL;
      } else {
        hi = pm + 1;
        break;
      }
    } else {
      hi = pm;
    }
  }

  if (ref->plen != lo) {
    ref->hash = seed;
    ref->plen = lo;
  }
  if (ret == NULL)
    ret = wormhole_hmap_get_kref(hmap, ref);
  debug_assert(ret);
  // ref->plen is the current depth, will be used in down()
  // ref now contains the prefix of meta-root
  return ret;
}

// pbuf: the current matched prefix of node
// klen0: the real key's len
  static struct wormleaf *
wormhole_meta_down(const struct wormhmap * const hmap, const struct wormkref * const ref, const struct wormmeta * const meta)
{
  struct wormleaf * ret;
  if (ref->plen < ref->klen) {
    const u32 id0 = ref->key[ref->plen];
    debug_assert(meta->bitmin != id0);
    if (meta->bitmin > id0) {
      ret = meta->lmost;
      if (meta->bitcount)
        ret = ret->prev;
    } else {
      const u32 id1 = wormhole_meta_bm_lt(meta, id0);
      const struct wormmeta * const child = wormhole_hmap_get_kref1(hmap, ref, id1);
      ret = child->rmost;
    }
  } else {
    // plen == klen
    debug_assert(ref->plen == ref->klen);
    ret = meta->lmost;
    if (ret->klen > ref->plen)
      ret = ret->prev;
  }

  // will lock
  __builtin_prefetch(&(ret->leaflock), 1, 1);
  return ret;
}

  static inline struct wormleaf *
wormhole_jump_leaf(const struct wormhmap * const hmap, const struct kv * const key)
{
  struct wormkref ref = {.hash = key->hash, .plen = key->klen, .klen = key->klen, .key = key->kv};

  const struct wormmeta * const meta = wormhole_meta_up(hmap, &ref);
  return wormhole_meta_down(hmap, &ref, meta);
}

  static inline struct wormleaf *
wormhole_jump_leaf_read(struct wormref * const ref, const struct kv * const key)
{
  struct wormhole * const map = ref->map;
  do {
    const struct wormhmap * const hmap = map->hmap;
    ref->qstate = (u64)hmap;
    struct wormleaf * const leaf = wormhole_jump_leaf(hmap, key);
    const u64 v = hmap->version;
    do {
      if (rwlock_trylock_read_nr(&(leaf->leaflock), 64)) {
        if (leaf->version <= v)
          return leaf;
        rwlock_unlock_read(&(leaf->leaflock));
      }
      ref->qstate = (u64)(map->hmap);
    } while (leaf->version <= v);
  } while (true);
}

  static inline struct wormleaf *
wormhole_jump_leaf_write(struct wormref * const ref, const struct kv * const key)
{
  struct wormhole * const map = ref->map;
  do {
    const struct wormhmap * const hmap = map->hmap;
    ref->qstate = (u64)hmap;
    struct wormleaf * const leaf = wormhole_jump_leaf(hmap, key);
    const u64 v = hmap->version;
    do {
      if (rwlock_trylock_write_nr(&(leaf->leaflock), 64)) {
        if (leaf->version <= v)
          return leaf;
        rwlock_unlock_write(&(leaf->leaflock));
      }
      ref->qstate = (u64)(map->hmap);
    } while (leaf->version <= v);
  } while (true);
}
// }}} jump

// leaf-only {{{
// bisect the first in eh[] whose hash is >= the given key's hash
  static u32
wormhole_leaf_bisect_hash(const struct wormleaf * const leaf, const u64 pkey)
{
  // bisect on the eh array (sorted by hash order)
  const u32 nr_keys = leaf->nr_keys;
  if (nr_keys == 0u)
    return 0u;

  u32 i = (((u32)pkey) * nr_keys) >> 16; // i < KPN
  debug_assert(i < nr_keys);
  const struct kv_entry * const eh = leaf->eh;
  if (pkey > eh[i].pkey) { // go right
    while ((i < nr_keys) && (pkey > eh[i].pkey))
      i++;
  } else { // go left
    while (i && (pkey <= eh[i - 1u].pkey))
      i--;
  }
  return i;
}

// find a matching key forward
  static u32
wormhole_leaf_seek_hash(const struct wormleaf * const leaf, const struct kv * const key,
    const u64 start, const u64 pkey)
{
  const u32 nr_keys = leaf->nr_keys;
  for (u32 i = start; i < nr_keys; i++) {
    if (leaf->eh[i].pkey != pkey)
      return nr_keys;

    if (kv_keymatch(key, u64_to_kv(leaf->eh[i].ptr)))
      return i;
  }
  return nr_keys;
}

// fast point-lookup
  static struct kv *
wormhole_leaf_match(const struct wormleaf * const leaf, const struct kv * const key)
{
  const u32 nr_keys = leaf->nr_keys;
  if (nr_keys == 0u)
    return NULL;

  const u64 pkey = kvmap_pkey(key->hash);
  u32 i = (pkey * nr_keys) >> 16; // i < KPN
  debug_assert(i < nr_keys);

  const struct kv_entry * const eh = leaf->eh;
  if (pkey > eh[i].pkey) { // go right, skip smaller
    while ((i < nr_keys) && (pkey > eh[i].pkey))
      i++;
  } else { // go left to first >=
    while (i && (pkey <= eh[i - 1u].pkey))
      i--;
  }
  while ((i < nr_keys) && (eh[i].pkey == pkey)) {
    struct kv * const curr = u64_to_kv(eh[i].ptr);
    if (kv_keymatch(key, curr))
      return curr;

    i++;
  }
  return NULL;
}
// }}} leaf-only

// get/probe {{{
  struct kv *
wormhole_get(struct wormref * const ref, const struct kv * const key, struct kv * const out)
{
  struct wormleaf * const leaf = wormhole_jump_leaf_read(ref, key);
  struct kv * tmp = wormhole_leaf_match(leaf, key);
  if (tmp) // found
    tmp = kv_dup2(tmp, out);
  rwlock_unlock_read(&(leaf->leaflock));
  return tmp;
}

  bool
wormhole_probe(struct wormref * const ref, const struct kv * const key)
{
  struct wormleaf * const leaf = wormhole_jump_leaf_read(ref, key);
  struct kv * tmp = wormhole_leaf_match(leaf, key);
  rwlock_unlock_read(&(leaf->leaflock));
  return tmp ? true : false;
}

  struct kv *
wormhole_get_unsafe(struct wormhole * const map, const struct kv * const key, struct kv * const out)
{
  struct wormleaf * const leaf = wormhole_jump_leaf(map->hmap, key);
  struct kv * tmp = wormhole_leaf_match(leaf, key);
  if (tmp) // found
    tmp = kv_dup2(tmp, out);
  return tmp;
}

  bool
wormhole_probe_unsafe(struct wormhole * const map, const struct kv * const key)
{
  struct wormleaf * const leaf = wormhole_jump_leaf(map->hmap, key);
  return wormhole_leaf_match(leaf, key) ? true : false;
}
// }}} get/probe

// }}} get

// set {{{

// leaf-only {{{
  static void
wormhole_leaf_sort_m2(struct kv_entry * const es, const u32 n1, const u32 n2)
{
  if (n1 == 0 || n2 == 0)
    return; // no need to sort

  struct kv_entry er[n1 + n2];
  struct kv_entry * ei = er;
  u32 r1 = n1;
  u32 r2 = n2;
  struct kv_entry * e1 = es;
  struct kv_entry * e2 = &(es[n1]);
  while (r1 && r2) {
    const int cmp = kv_keycompare(u64_to_kv(e1->ptr), u64_to_kv(e2->ptr));
    if (cmp < 0) {
      *ei = *e1; e1++; r1--;
    } else if (cmp > 0) {
      *ei = *e2; e2++; r2--;
    } else {
      debug_die();
    }
    ei++;
  }
  const u32 done = n1 + n2 - r1 - r2;
  if (r1)
    memmove(es + done, e1, sizeof(*es) * r1);

  memcpy(es, er, sizeof(*es) * done);
}

// make sure all keys are sorted in a leaf node
  static void
wormhole_leaf_sync_sorted(struct wormleaf * const leaf)
{
  const u32 s = leaf->nr_sorted;
  const u32 n = leaf->nr_keys;
  if (s == n)
    return;

  if ((n < 8) || (s < 4u) || (((s - 2u) << 2) < n)) { // too few sorted
    __kvmap_entry_qsort(leaf->es, n);
  } else { // worth a two-step sort
    __kvmap_entry_qsort(&(leaf->es[s]), n - s);
    // merge-sort inplace
    wormhole_leaf_sort_m2(leaf->es, s, (n - s));
  }
  leaf->nr_sorted = n;
}

// test if leaf->es[ip] is a prefix of leaf->es[ik]
  static bool
wormhole_split_can_cut(const struct kv_entry * const es, const u32 i1, const u32 i2)
{
  debug_assert(i1 < i2);
  const struct kv * const k1 = u64_to_kv(es[i1].ptr);
  const struct kv * const k2 = u64_to_kv(es[i2].ptr);
  // rule of no-cut: k1 and k2 have the same TTK.
  // first test if k1 is a prefix of k2, then trim k to see if it has a zero-tail.
  if (kv_key_is_prefix(k1, k2)) {
    u32 tklen = k2->klen;
    while ((tklen > k1->klen) && (k2->kv[tklen - 1u] == 0u)) tklen--;
    if (tklen <= k1->klen)
      return false;
  }
  return true;
}

// assumes i1/i2 has been verified by split_can_cut
  static u32
wormhole_split_cut_alen(const struct kv_entry * const es, const u32 i1, const u32 i2)
{
  debug_assert(i1 < i2);
  struct kv * const k1 = u64_to_kv(es[i1].ptr);
  struct kv * const k2 = u64_to_kv(es[i2].ptr);
  u32 alen = kv_key_lcp(k1, k2) + 1u;
  while ((alen < k2->klen) && (k2->kv[alen - 1u] == 0u)) alen++;
  debug_assert(k2->kv[alen - 1u]);
  return alen;
}

  static bool
wormhole_split_can_cut_alen(const struct kv_entry * const es, const u32 i1, const u32 i2, const u32 alen)
{
  debug_assert(i1 < i2);
  struct kv * const k1 = u64_to_kv(es[i1].ptr);
  struct kv * const k2 = u64_to_kv(es[i2].ptr);
  const u8 c1 = (k1->klen < alen) ? 0u : k1->kv[alen - 1u];
  const u8 c2 = (k2->klen < alen) ? 0u : k2->kv[alen - 1u];
  return (c1 == c2) ? false : true;
}

// determine where to cut at leaf
  static u32
wormhole_split_cut(const struct wormleaf * const leaf)
{
  debug_assert(leaf->nr_keys == WORMHOLE_KPN);
  debug_assert(leaf->nr_sorted == WORMHOLE_KPN);
  const struct kv_entry * const es = leaf->es;
  u32 lo = 0;
  u32 hi = WORMHOLE_KPN - 1u;
  if (wormhole_split_can_cut(es, lo, hi) == false)
    return WORMHOLE_KPN;

  const u32 alen = wormhole_split_cut_alen(es, lo, hi);
  while ((lo + 1u) < hi) {
    const u32 mid = (lo + hi + 1u) >> 1u;
    if (mid <= WORMHOLE_MID) { // try right
      if (wormhole_split_can_cut_alen(es, mid, hi, alen))
        lo = mid;
      else
        hi = mid;
    } else { // try left
      if (wormhole_split_can_cut_alen(es, lo, mid, alen))
        hi = mid;
      else
        lo = mid;
    }
  }
  return hi;
}

  static inline void
wormhole_leaf_sync_hash(struct wormleaf * const leaf)
{
  // sort eh in hash order
  qsort(leaf->eh, leaf->nr_keys, sizeof(leaf->eh[0]), __kvmap_entry_pkey_comp);
}

  static void
wormhole_split_leaf_move(struct wormleaf * const leaf1, struct wormleaf * const leaf2, const u32 cut)
{
  const u32 nr_all = leaf1->nr_keys;
  const size_t s1 = sizeof(leaf1->es[0]) * cut;
  const size_t s2 = sizeof(leaf2->es[0]) * (nr_all - cut);
  // move es
  memcpy(leaf2->es, &(leaf1->es[cut]), s2);
  // duplicate to eh
  memcpy(leaf1->eh, leaf1->es, s1);
  memcpy(leaf2->eh, leaf2->es, s2);
  // metadata
  leaf1->nr_keys = cut;
  leaf1->nr_sorted = cut;
  leaf2->nr_keys = (nr_all - cut);
  leaf2->nr_sorted = (nr_all - cut);
  wormhole_leaf_sync_hash(leaf1);
  wormhole_leaf_sync_hash(leaf2);
}

  static void
wormhole_set_update(struct wormhole * const map, struct wormleaf * const leaf, const u32 im, const struct kv * const new)
{
  // locate in es (is)
  const u64 v64 = leaf->eh[im].v64;
  const u32 nr = leaf->nr_keys;
  u32 is;
  for (is = 0; is < nr; is++) {
    if (leaf->es[is].v64 == v64)
      break;
  }
  debug_assert(is < nr);

  kvmap_put_entry(&(map->mm), &(leaf->eh[im]), new);
  leaf->es[is] = leaf->eh[im];
}

  static void
wormhole_set_insert(struct wormleaf * const leaf, const u32 ii, const struct kv * const new, const u64 pkey)
{
  const u32 nr0 = leaf->nr_keys;
  // insert into leaf
  // append to es (delayed sort)
  leaf->es[nr0].ptr = kv_to_u64(new);
  leaf->es[nr0].pkey = pkey;

  // insert into eh (use ii)
  memmove(&(leaf->eh[ii + 1]), &(leaf->eh[ii]), sizeof(leaf->eh[ii]) * (nr0 - ii));
  leaf->eh[ii] = leaf->es[nr0];
  leaf->nr_keys++;
}

// create an anchor for leaf2
  static void
wormhole_split_anchor(struct wormleaf * const leaf1, struct wormleaf * const leaf2)
{
  struct kv * const key1 = u64_to_kv(leaf1->es[leaf1->nr_sorted - 1lu].ptr);
  struct kv * const key2 = u64_to_kv(leaf2->es[0].ptr);
  const u32 key2len = key2->klen;
  u32 alen = kv_key_lcp(key1, key2) + 1lu;

  // anchor must end with non-zero
  while ((alen < key2len) && (key2->kv[alen - 1u] == 0u)) alen++;
  debug_assert(alen <= key2len);

  // now we have the correct alen
  struct kv * const anchor2 = wormhole_alloc_key(alen);
  debug_assert(anchor2);
  kv_refill(anchor2, key2->kv, alen, NULL, 0);
  leaf2->anchor = anchor2;
  leaf2->klen = anchor2->klen;
}
// }}} leaf-only

// split {{{
// move keys starting with [cut] in leaf1 to leaf2
  static struct wormleaf *
wormhole_split_leaf(struct wormhole * const map, struct wormleaf * const leaf1, const u32 cut)
{
  // create leaf (2)
  struct wormleaf * const leaf2 = wormhole_alloc_leaf(map, leaf1, leaf1->next, NULL);
  wormhole_split_leaf_move(leaf1, leaf2, cut);

  // anchor of leaf2
  wormhole_split_anchor(leaf1, leaf2);
  rwlock_lock_write(&(leaf2->leaflock));
  return leaf2;
}

  static void
wormhole_split_meta_new(struct wormhmap * const hmap, const struct kv * const pbuf,
    struct wormleaf * const leaf, const bool bitmap)
{
  // create a new node
  struct kv * const anchor = wormhole_alloc_key(pbuf->klen);
  kv_dup2_key(pbuf, anchor);
  struct wormmeta * const meta = wormhole_alloc_meta(hmap->map, leaf, anchor);

  if (bitmap)
    wormhole_meta_bm_set(meta, pbuf->kv[pbuf->klen]);

  wormhole_hmap_set(hmap, anchor, meta);
}

  static void
wormhole_split_meta_extend(struct wormhmap * const hmap, struct wormmeta * const meta,
    const struct kv * const pbuf)
{
  debug_assert(meta->lmost == meta->rmost);
  wormhole_meta_bm_set(meta, 0u);
  // new anchor at +1
  struct kv * const anchor1 = wormhole_alloc_key(pbuf->klen + 1u);
  kv_dup2_key(pbuf, anchor1);
  anchor1->kv[anchor1->klen] = 0u;
  wormhole_prefix_inc_short(anchor1, anchor1->klen + 1u);
  struct wormmeta * const meta1 = wormhole_alloc_meta(hmap->map, meta->lmost, anchor1);
  wormhole_hmap_set(hmap, anchor1, meta1);
}

// return true if a new node is created
  static bool
wormhole_split_meta_one(struct wormhmap * const hmap, const struct kv * const pbuf,
    struct wormleaf * const leaf, const bool setchild)
{
  struct wormmeta * const meta = wormhole_hmap_get(hmap, pbuf);
  if (meta == NULL) {
    wormhole_split_meta_new(hmap, pbuf, leaf, setchild);
    return true;
  }

  // push down leaf
  if (meta->bitcount == 0u) {
    wormhole_split_meta_extend(hmap, meta, pbuf);
  }

  // mark leaf's child bit
  wormhole_meta_bm_set(meta, pbuf->kv[pbuf->klen]);

  // lmost rmost
  if (meta->lmost == leaf->next)
    meta->lmost = leaf;

  if (meta->rmost == leaf->prev)
    meta->rmost = leaf;
  return false;
}

// for leaf1, a leaf2 is already linked at its right side.
// this function updates the meta-map by moving leaf1 and hooking leaf2 at correct positions
  static void
wormhole_split_meta_hmap(struct wormhmap * const hmap, struct wormleaf * const leaf)
{
  // pbuf
  const struct kv * const anchor = leaf->anchor;
  const u32 len2 = leaf->next ? leaf->next->anchor->klen : 0lu;
  const u32 maxplen = anchor->klen > len2 ? anchor->klen : len2;
  struct kv * const pbuf = alloca(sizeof(struct kv) + maxplen + 1);
  kv_dup2_key(anchor, pbuf);
  // child bit will be set base on this value

  // left branches
  const u32 lcp1 = kv_key_lcp(leaf->prev->anchor, leaf->anchor);
  const u32 lcp2 = leaf->next ? kv_key_lcp(leaf->anchor, leaf->next->anchor) : 0u;
  u32 i = (lcp1 < lcp2) ? lcp1 : lcp2;

  wormhole_prefix(pbuf, i);
  while (i < anchor->klen) {
    wormhole_split_meta_one(hmap, pbuf, leaf, true);
    i++;
    wormhole_prefix_inc_short(pbuf, i);
  }
  while (true) {
    debug_assert(i <= maxplen);
    pbuf->kv[i] = 0;
    if (wormhole_split_meta_one(hmap, pbuf, leaf, false)) break;
    i++;
    wormhole_prefix_inc_short(pbuf, i);
  }

  // adjust maxplen
  if (i > hmap->maxplen)
    hmap->maxplen = i;
}

  static inline struct wormhmap *
wormhole_hmap_sibling(struct wormhole * const map, struct wormhmap * const hmap0)
{
  return (&(map->hmap2[0])) == hmap0 ? (&(map->hmap2[1])) : (&(map->hmap2[0]));
}

// all locks will be released before returning
  static void
wormhole_split_meta_ref(struct wormref * const ref, struct wormleaf * const leaf2)
{
  struct wormhole * const map = ref->map;
  // metalock
  while (rwlock_trylock_write_nr(&(map->metalock), 64) == false)
    ref->qstate = (u64)(map->hmap);

  struct wormhmap * const hmap0 = map->hmap;
  struct wormhmap * const hmap1 = wormhole_hmap_sibling(map, hmap0);
  const u64 v1 = hmap0->version + 1;

  // new versions
  struct wormleaf * const leaf1 = leaf2->prev;

  // link
  leaf1->next = leaf2;
  if (leaf2->next)
    leaf2->next->prev = leaf2;
  // update versions
  leaf1->version = v1;
  leaf2->version = v1;
  hmap1->version = v1;

  wormhole_split_meta_hmap(hmap1, leaf2);
  ref->qstate = (u64)(hmap1);
  // switch hmap
  cpu_cfence();
  map->hmap = hmap1;
  cpu_cfence();

  rwlock_unlock_write(&(leaf1->leaflock));
  rwlock_unlock_write(&(leaf2->leaflock));

  qsbr_wait(map->qsbr, (u64)hmap1);

  wormhole_split_meta_hmap(hmap0, leaf2);
  rwlock_unlock_write(&(map->metalock));
}

// all locks will be released before returning
// leaf1->lock is already taken
  static bool
wormhole_split_insert(struct wormref * const ref, struct wormleaf * const leaf1,
    const struct kv * const new, const u64 pkey)
{
  wormhole_leaf_sync_sorted(leaf1);
  // corner case that we don't handle for now.
  // TODO: we can still split but don't change any metadata.
  //       The split right leaf on the chain can be a shadow companion for the left leaf.
  //       But that will make the updates to lmost/rmost more complex.
  const u32 cut = wormhole_split_cut(leaf1);
  debug_assert(cut && (cut < leaf1->nr_keys));

  struct wormhole * const map = ref->map;
  struct wormleaf * const leaf2 = wormhole_split_leaf(map, leaf1, cut);

  // the compare must regard the anchor's extensions
  const int cmp = kv_keycompare(new, leaf2->anchor);
  struct wormleaf * const leaf = (cmp < 0) ? leaf1 : leaf2;
  const u32 ii = wormhole_leaf_bisect_hash(leaf, pkey);
  wormhole_set_insert(leaf, ii, new, pkey);

  wormhole_split_meta_ref(ref, leaf2);
  return true;
}

  static void
wormhole_split_meta_unsafe(struct wormhole * const map, struct wormleaf * const leaf2)
{
  struct wormhmap * const hmap0 = map->hmap;
  struct wormhmap * const hmap1 = wormhole_hmap_sibling(map, hmap0);

  // new versions
  struct wormleaf * const leaf1 = leaf2->prev;

  // link
  leaf1->next = leaf2;
  if (leaf2->next)
    leaf2->next->prev = leaf2;

  wormhole_split_meta_hmap(hmap1, leaf2);
  wormhole_split_meta_hmap(hmap0, leaf2);
}

  static bool
wormhole_split_insert_unsafe(struct wormhole * const map, struct wormleaf * const leaf1,
    const struct kv * const new, const u64 pkey)
{
  wormhole_leaf_sync_sorted(leaf1);
  // corner case that we don't handle for now.
  // TODO: we can still split but don't change any metadata.
  //       The split right leaf on the chain can be a shadow companion for the left leaf.
  //       But that will make the updates to lmost/rmost more complex.
  const u32 cut = wormhole_split_cut(leaf1);
  debug_assert(cut && (cut < leaf1->nr_keys));

  struct wormleaf * const leaf2 = wormhole_split_leaf(map, leaf1, cut);

  // the compare must regard the anchor's extensions
  const int cmp = kv_keycompare(new, leaf2->anchor);
  struct wormleaf * const leaf = (cmp < 0) ? leaf1 : leaf2;
  const u32 ii = wormhole_leaf_bisect_hash(leaf, pkey);
  wormhole_set_insert(leaf, ii, new, pkey);

  wormhole_split_meta_unsafe(map, leaf2);
  return true;
}
// }}} split

// set {{{
  bool
wormhole_set(struct wormref * const ref, const struct kv * const kv0)
{
  struct wormhole * const map = ref->map;
  struct kv * const new = wormhole_alloc_kv(map, kv0->klen, kv0->vlen);
  if (new == NULL)
    return false;
  kv_dup2(kv0, new);

  struct wormleaf * const leaf = wormhole_jump_leaf_write(ref, new);
  // insertion point
  const u64 pkey = kvmap_pkey(new->hash);
  const u32 ii = wormhole_leaf_bisect_hash(leaf, pkey);
  // maybe a match
  const u32 im = wormhole_leaf_seek_hash(leaf, new, ii, pkey);
  if (im < leaf->nr_keys) { // overwrite
    wormhole_set_update(map, leaf, im, new);
    rwlock_unlock_write(&(leaf->leaflock));
    return true;
  }
  // insert
  if (leaf->nr_keys < WORMHOLE_KPN) { // just insert
    wormhole_set_insert(leaf, ii, new, pkey);
    rwlock_unlock_write(&(leaf->leaflock));
    return true;
  }

  // changes hmap
  return wormhole_split_insert(ref, leaf, new, pkey);
}

  bool
wormhole_set_unsafe(struct wormhole * const map, const struct kv * const kv0)
{
  struct kv * const new = wormhole_alloc_kv(map, kv0->klen, kv0->vlen);
  if (new == NULL)
    return false;
  kv_dup2(kv0, new);

  struct wormleaf * const leaf = wormhole_jump_leaf(map->hmap, new);
  // insertion point
  const u64 pkey = kvmap_pkey(new->hash);
  const u32 ii = wormhole_leaf_bisect_hash(leaf, pkey);
  // maybe a match
  const u32 im = wormhole_leaf_seek_hash(leaf, new, ii, pkey);
  if (im < leaf->nr_keys) { // overwrite
    wormhole_set_update(map, leaf, im, new);
    return true;
  }
  // insert
  if (leaf->nr_keys < WORMHOLE_KPN) { // just insert
    wormhole_set_insert(leaf, ii, new, pkey);
    return true;
  }

  // changes hmap
  return wormhole_split_insert_unsafe(map, leaf, new, pkey);
}
// }}} set

// }}} set

// del {{{
  static void
wormhole_leaf_del_one(struct wormhole * const map, struct wormleaf * const leaf, const u32 ih)
{
  const u32 nr_keys = leaf->nr_keys;
  // remove it from sorted
  u32 is;
  for (is = 0; is < nr_keys; is++) {
    if (leaf->es[is].v64 == leaf->eh[ih].v64) {
      if (is < (nr_keys - 1u))
        leaf->es[is] = leaf->es[nr_keys - 1u];

      break;
    }
  }
  debug_assert(is < nr_keys);
  if (leaf->nr_sorted > is)
    leaf->nr_sorted = is;

  kvmap_put_entry(&(map->mm), &(leaf->eh[ih]), NULL);
  // eh: shift left
  memmove(&(leaf->eh[ih]), &(leaf->eh[ih + 1u]), sizeof(leaf->eh[ih]) * (nr_keys - ih - 1u));
  leaf->nr_keys--;
}

// all go to leaf1
  static void
wormhole_merge_leaf_move(struct wormleaf * const leaf1, struct wormleaf * const leaf2)
{
  const u32 nr1 = leaf1->nr_keys;
  const u32 nr2 = leaf2->nr_keys;
  debug_assert((nr1 + nr2) <= WORMHOLE_KPN);
  memcpy(&(leaf1->eh[nr1]), &(leaf2->eh[0]), sizeof(leaf2->eh[0]) * nr2);
  memcpy(&(leaf1->es[nr1]), &(leaf2->eh[0]), sizeof(leaf2->eh[0]) * nr2);
  leaf1->nr_keys = nr1 + nr2; // nr_sorted remain unchanged
  if (nr1 && nr2)
    wormhole_leaf_sync_hash(leaf1);
}

  static void
wormhole_merge_meta(struct wormhmap * const hmap, struct wormleaf * const leaf)
{
  struct wormleaf * const prev = leaf->prev;
  struct wormleaf * const next = leaf->next;
  const struct kv * const anchor0 = leaf->anchor;
  const u32 lcp1 = prev ? kv_key_lcp(prev->anchor, anchor0) : 0;
  const u32 lcp2 = next ? kv_key_lcp(next->anchor, anchor0) : 0;
  const u32 maxplen = hmap->maxplen;
  const u64 bsize = sizeof(struct kv) + maxplen;
  struct kv * const pbuf = alloca(bsize);
  kv_dup2_key(anchor0, pbuf);
  u32 i = lcp1 < lcp2 ? lcp1 : lcp2;
  // lmost & rmost
  struct wormmeta * parent = NULL;
  wormhole_prefix(pbuf, i);
  while (true) {
    debug_assert(i <= maxplen);
    struct wormmeta * const meta = wormhole_hmap_get(hmap, pbuf);
    debug_assert(meta);
    if (meta->lmost == meta->rmost) {
      debug_assert(meta->lmost == leaf);
      const u32 cnt = meta->bitcount;
      wormhole_hmap_del(hmap, pbuf);
      wormhole_free_key(meta->key);
      slab_free(hmap->map->slab_meta, meta);
      if (parent) {
        wormhole_meta_bm_clear(parent, pbuf->kv[i - 1u]);
        parent = NULL;
      }
      if (cnt == 0u) break;
    } else {
      if (meta->lmost == leaf)
        meta->lmost = next;

      if (meta->rmost == leaf)
        meta->rmost = prev;
      parent = meta;
    }
    if (i >= anchor0->klen)
      pbuf->kv[i] = 0u;
    i++;
    wormhole_prefix_inc_short(pbuf, i);
  }
}

// all locks will be released before returning
// merge leaf2 to leaf1, removing all metadata to leaf2 and leaf2 itself
  static void
wormhole_merge_ref(struct wormref * const ref, struct wormleaf * const leaf1, struct wormleaf * const leaf2)
{
  wormhole_merge_leaf_move(leaf1, leaf2);

  struct wormhole * const map = ref->map;
  struct wormhmap * const hmap0 = map->hmap;
  struct wormhmap * const hmap1 = wormhole_hmap_sibling(map, hmap0);
  const u64 v1 = hmap0->version + 1;

  leaf1->next = leaf2->next;
  if (leaf2->next)
    leaf2->next->prev = leaf1;

  leaf1->version = v1;
  leaf2->version = v1;
  hmap1->version = v1;

  wormhole_merge_meta(hmap1, leaf2);
  ref->qstate = (u64)(hmap1);

  cpu_cfence();
  map->hmap = hmap1;
  cpu_cfence();

  rwlock_unlock_write(&(leaf1->leaflock));
  rwlock_unlock_write(&(leaf2->leaflock));

  qsbr_wait(map->qsbr, (u64)hmap1);

  wormhole_merge_meta(hmap0, leaf2);
  rwlock_unlock_write(&(map->metalock));
  // leaf2 is now safe to be removed
  wormhole_free_key(leaf2->anchor);
  slab_free(map->slab_leaf, leaf2);
}

// all locks will be released before returning
// merge may fail if (1) too large combined; (2) didn't acquire lock before timeout
  static bool
wormhole_try_merge(struct wormref * const ref, struct wormleaf * const leaf)
{
  struct wormhole * const map = ref->map;
  while (rwlock_trylock_write_nr(&(map->metalock), 64) == false)
    ref->qstate = (u64)(map->hmap);

  // now leaf and meta are both locked (w)
  struct wormleaf * const next = leaf->next;
  debug_assert(next);
  if (rwlock_trylock_write_nr(&(next->leaflock), 64)) {
    // three locked
    if ((leaf->nr_keys + next->nr_keys) <= WORMHOLE_KPN) {
      wormhole_merge_ref(ref, leaf, next);
      return true;
    }
    rwlock_unlock_write(&(next->leaflock));
  }

  rwlock_unlock_write(&(map->metalock));
  rwlock_unlock_write(&(leaf->leaflock));
  return false;
}

  static bool
wormhole_force_merge(struct wormref * const ref, struct wormleaf * const leaf)
{
  struct wormhole * const map = ref->map;
  while (rwlock_trylock_write_nr(&(map->metalock), 64) == false)
    ref->qstate = (u64)(map->hmap);

  // now leaf and meta are both locked (w)
  struct wormleaf * const next = leaf->next;
  debug_assert(next);
  while (rwlock_trylock_write_nr(&(next->leaflock), 64) == false)
    ref->qstate = (u64)(map->hmap);
  // three locked
  wormhole_merge_ref(ref, leaf, next);
  return true;
}

  bool
wormhole_del(struct wormref * const ref, const struct kv * const key)
{
  struct wormleaf * const leaf = wormhole_jump_leaf_write(ref, key);
  const u64 pkey = kvmap_pkey(key->hash);
  const u32 ii = wormhole_leaf_bisect_hash(leaf, pkey);
  const u32 im = wormhole_leaf_seek_hash(leaf, key, ii, pkey);
  bool r = false;
  if (im < leaf->nr_keys) { // found
    wormhole_leaf_del_one(ref->map, leaf, im);
    r = true;

    // maybe merge
    if ((leaf->nr_keys == 0u) && leaf->next) {
      (void)wormhole_force_merge(ref, leaf);
      return r;
    }
    struct wormleaf * const next = leaf->next;
    if (next && ((leaf->nr_keys + next->nr_keys) < WORMHOLE_KPN_MRG)) {
      (void)wormhole_try_merge(ref, leaf);
      return r;
    }
  }

  rwlock_unlock_write(&(leaf->leaflock));
  return r;
}

  static void
wormhole_merge_unsafe(struct wormhole * const map, struct wormleaf * const leaf1, struct wormleaf * const leaf2)
{
  wormhole_merge_leaf_move(leaf1, leaf2);
  struct wormhmap * const hmap0 = map->hmap;
  struct wormhmap * const hmap1 = wormhole_hmap_sibling(map, hmap0);

  leaf1->next = leaf2->next;
  if (leaf2->next)
    leaf2->next->prev = leaf1;
  wormhole_merge_meta(hmap1, leaf2);
  wormhole_merge_meta(hmap0, leaf2);
  wormhole_free_key(leaf2->anchor);
  slab_free(map->slab_leaf, leaf2);
}

  bool
wormhole_del_unsafe(struct wormhole * const map, const struct kv * const key)
{
  struct wormleaf * const leaf = wormhole_jump_leaf(map->hmap, key);
  const u64 pkey = kvmap_pkey(key->hash);
  const u32 ii = wormhole_leaf_bisect_hash(leaf, pkey);
  const u32 im = wormhole_leaf_seek_hash(leaf, key, ii, pkey);
  if (im < leaf->nr_keys) { // found
    wormhole_leaf_del_one(map, leaf, im);

    struct wormleaf * const next = leaf->next;
    if (next && ((leaf->nr_keys + next->nr_keys) < WORMHOLE_KPN_MRG)) {
      wormhole_merge_unsafe(map, leaf, next);
    }
    return true;
  }
  return false;
}
// }}} del

// misc {{{
  inline struct wormref *
wormhole_ref(struct wormhole * const map)
{
  struct wormref * const ref = yalloc(sizeof(*ref));
  ref->qstate = 0lu;
  ref->map = map;
  const bool r = qsbr_register(map->qsbr, &(ref->qstate));
  (void)r;
  debug_assert(r);
  return ref;
}

  inline struct wormhole *
wormhole_unref(struct wormref * const ref)
{
  struct wormhole * const map = ref->map;
  qsbr_unregister(map->qsbr, &(ref->qstate));
  free(ref);
  return map;
}

  static void
wormhole_clean1(struct wormhole * const map)
{
  // meta
  for (u64 x = 0; x < 2; x++) {
    const u64 nr_slots = map->hmap2[x].mask + 1lu;
    for (u64 s = 0; s < nr_slots; s++) {
      struct kvmap_slot * const slot = &(map->hmap2[x].pmap[s]);
      for (u64 i = 0; i < KVMAP_SLOT_NR; i++) {
        struct kv_entry * const e = &(slot->e[i]);
        if (e->v64 == 0lu)
          continue;
        struct wormmeta * const meta = __u64_to_meta(e->ptr);
        wormhole_free_key(meta->key);
        slab_free(map->slab_meta, meta);
        e->v64 = 0lu;
        map->hmap2[x].wmap[s].t[i] = 0u;
      }
    }
    map->hmap2[x].maxplen = 0u;
  }
  // leaf
  struct wormleaf * leaf = map->leaf0;
  while (leaf) {
    struct wormleaf * const next = leaf->next;
    wormhole_free_key(leaf->anchor);
    for (u64 i = 0; i < leaf->nr_keys; i++)
      kvmap_put_entry(&(map->mm), &(leaf->eh[i]), NULL);
    slab_free(map->slab_leaf, leaf);
    leaf = next;
  }
  map->leaf0 = NULL;
}

  void
wormhole_clean(struct wormhole * const map)
{
  wormhole_clean1(map);
  wormhole_create_leaf0(map);
}

  void
wormhole_destroy(struct wormhole * const map)
{
  //wormhole_verify(map);
  wormhole_clean1(map);
  for (u64 x = 0; x < 2; x++) {
    pages_unmap(map->hmap2[x].pmap, map->hmap2[x].psize);
    pages_unmap(map->hmap2[x].wmap, map->hmap2[x].wsize);
  }
  qsbr_destroy(map->qsbr);
  slab_destroy(map->slab_meta);
  slab_destroy(map->slab_leaf);
  free(map);
}

  struct wormhole_iter *
wormhole_iter_create(struct wormref * const ref)
{
  struct wormhole_iter * const iter = malloc(sizeof(*iter));
  iter->ref = ref;
  iter->leaf = NULL;
  iter->next_id = 0;
  wormhole_iter_seek(iter, NULL);
  return iter;
}

// assumes there in no duplicated keys
// bisect the first key that is >= the given key
  static u32
wormhole_leaf_bisect_sorted(const struct wormleaf * const leaf, const struct kv * const key)
{
  u32 lo = 0u;
  u32 hi = leaf->nr_sorted;
  while (lo < hi) {
    u32 i = (lo + hi) >> 1;
    const int cmp = kv_keycompare(u64_to_kv(leaf->es[i].ptr), key);
    if (cmp < 0)  //  [i] < key
      lo = i + 1u;
    else if (cmp > 0)
      hi = i;
    else // same key
      return i;
  }
  return lo;
}

  void
wormhole_iter_seek(struct wormhole_iter * const iter, const struct kv * const key)
{
  if (iter->leaf)
    rwlock_unlock_read(&(iter->leaf->leaflock));
  struct wormref * const ref = iter->ref;
  struct wormhole * const map = ref->map;

  if (key == NULL) {
    struct wormleaf * const leaf0 = iter->ref->map->leaf0;
    iter->leaf = leaf0;
    iter->next_id = 0;
    while (rwlock_trylock_read_nr(&(leaf0->leaflock), 64) == false)
      ref->qstate = (u64)(map->hmap);
    wormhole_leaf_sync_sorted(leaf0);
    return;
  }

  struct wormleaf * const leaf = wormhole_jump_leaf_read(ref, key);
  wormhole_leaf_sync_sorted(leaf);
  const u32 id = wormhole_leaf_bisect_sorted(leaf, key);
  if (id < leaf->nr_sorted) {
    iter->leaf = leaf;
    iter->next_id = id;
  } else {
    struct wormleaf * const next = leaf->next;
    iter->leaf = next;
    iter->next_id = 0;
    if (next) {
      while (rwlock_trylock_read_nr(&(next->leaflock), 64) == false)
        ref->qstate = (u64)(map->hmap);
      wormhole_leaf_sync_sorted(next);
    }
    rwlock_unlock_read(&(leaf->leaflock));
  }
}

  struct kv *
wormhole_iter_next(struct wormhole_iter * const iter, struct kv * const out)
{
  if (iter->leaf == NULL) return NULL;
  while (iter->next_id >= iter->leaf->nr_sorted) {
    struct wormleaf * const next = iter->leaf->next;
    if (next) {
      struct wormref * const ref = iter->ref;
      struct wormhole * const map = ref->map;
      while (rwlock_trylock_read_nr(&(next->leaflock), 64) == false)
        ref->qstate = (u64)(map->hmap);
    }
    rwlock_unlock_read(&(iter->leaf->leaflock));
    iter->leaf = next;
    if (iter->leaf == NULL) return NULL;
    iter->next_id = 0;
    wormhole_leaf_sync_sorted(iter->leaf);
  }
  debug_assert(iter->leaf);
  debug_assert(iter->next_id < iter->leaf->nr_sorted);
  struct kv * const kv = u64_to_kv(iter->leaf->es[iter->next_id].ptr);
  struct kv * const ret = kv_dup2(kv, out);
  iter->next_id++;
  return ret;
}

  void
wormhole_iter_destroy(struct wormhole_iter * const iter)
{
  if (iter->leaf)
    rwlock_unlock_read(&(iter->leaf->leaflock));
  free(iter);
}

// }}} misc

// }}} wormhole

// fdm: marker
