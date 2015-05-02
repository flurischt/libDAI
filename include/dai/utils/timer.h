#ifndef __defined_libdai_utils_timer_h
#define __defined_libdai_utils_timer_h

#include <cstdint>
#include <string>
#include <vector>
#include <iostream>
#include <sys/time.h>

/*
 * Sample code from lecture FNC
 * TODO: Give More credits
 */

/* ================ GNU C and possibly other UNIX compilers ================ */
#ifndef WIN32
   #if defined(__GNUC__) || defined(__linux__)
      #define VOLATILE __volatile__
      #define ASM __asm__
   #else
      /* if we're neither compiling with gcc or under linux, we can hope
       * the following lines work, they probably won't */
      #define ASM asm
      #define VOLATILE
   #endif

/* ================================= WIN32 ================================= */
#else

#endif

/* This is the RDTSC timer.
 * RDTSC is an instruction on several Intel and compatible CPUs that reads the
 * Time Stamp Counter. The Intel manuals contain more information.
 */
#define COUNTER_LO(a) ((a).int32.lo)
#define COUNTER_HI(a) ((a).int32.hi)
#define COUNTER_VAL(a) ((a).int64)

#define COUNTER(a) \
   ((unsigned long long)COUNTER_VAL(a))

#define COUNTER_DIFF(a,b) \
   (COUNTER(a)-COUNTER(b))

/* ==================== GNU C and possibly other UNIX compilers ===================== */
#ifndef WIN32

   typedef union
   {
      int64_t int64;
      struct {int32_t lo, hi;} int32;
   } tsc_counter;

  #define RDTSC(cpu_c) \
     ASM VOLATILE ("rdtsc" : "=a" ((cpu_c).int32.lo), "=d"((cpu_c).int32.hi))
   #define CPUID() \
      ASM VOLATILE ("cpuid" : : "a" (0) : "bx", "cx", "dx" )

/* ================================= WIN32 ================================= */
#else

   typedef union
   {
      int64_t int64;
      struct {int32_t lo, hi;} int32;
   } tsc_counter;

   #define RDTSC(cpu_c)  \
   {  \
      __asm rdtsc  \
      __asm mov (cpu_c).int32.lo,eax  \
      __asm mov (cpu_c).int32.hi,edx  \
   }

   #define CPUID() \
   { \
      __asm mov eax, 0 \
      __asm cpuid \
   }

#endif


class TimerAbstract {
public:
   virtual void tic() = 0;
   virtual double toc() = 0;
};

/*****************************************************************************/
class TimerGOD : public TimerAbstract{
public:
   TimerGOD() : TimerAbstract()
   {
      start_time.tv_sec  = 0;
      start_time.tv_usec = 0;
      stop_time.tv_sec   = 0;
      stop_time.tv_usec  = 0;
   }

   virtual void tic() {
      gettimeofday(&start_time, NULL);
   }

   virtual double toc() {
      gettimeofday(&stop_time, NULL);
      return (stop_time.tv_sec - start_time.tv_sec)
           + (stop_time.tv_usec - start_time.tv_usec)*1e-6;
   }

private:
   struct timeval start_time, stop_time;
};

/*****************************************************************************/
class TimerTSC : public TimerAbstract {
public:
   TimerTSC() : TimerAbstract()
   {
   }

   virtual void tic() {
      CPUID();
      RDTSC(start);
   }

   virtual double toc() {
      tsc_counter end;
      RDTSC(end);
      CPUID();
      return COUNTER_VAL(end) - COUNTER_VAL(start);
   }

private:
   tsc_counter start;
};

typedef TimerTSC Timer;

/*****************************************************************************/
struct PerformanceStats
{
   void addMeasurement(double duration,
                       const std::string &label,
                       std::size_t iterations = 1)
   {
      m_durations.push_back(duration);
      m_iterations.push_back(iterations);
      m_labels.push_back(label);
   }

   void print() const
   {
      const std::size_t n = m_durations.size();
      printf("------------------------------------------------\n");
      printf("Labels    :   Loops      Total    Average\n");
      for (std::size_t i=0; i<n; ++i)
      {
         printf("%-10s: %7zu %10g %10g\n", m_labels[i].c_str(),
                m_iterations[i], m_durations[i],
                m_durations[i]/m_iterations[i]);
      }
      printf("------------------------------------------------\n");
   }

   void clear()
   {
      m_durations.clear();
      m_iterations.clear();
      m_labels.clear();
   }

private:
   std::vector<double> m_durations;
   std::vector<std::size_t> m_iterations;
   std::vector<std::string> m_labels;
};

#endif //__defined_libdai_utils_timer_h
