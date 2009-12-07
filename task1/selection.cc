/*

  Selection algorithms in C++
  Copyright: Jyrki Katajainen 1999, 2009
  E-mail: jyrki@diku.dk

  The routines
    select_with_merge_sort,
    select_with_STL_sort,
    randomized_select,
    select, and
    nth_element from the STL
  do not only select the k'th largest of the given n elements
  but also partition the input into three parts: the elements
  smaller than or equal to the selected element, the selected
  element at the k'th location of the sequence, and the
  elements larger than or equal to the selected element.

  The item abstraction of the LEDA is used. According to
  the standard STL terminology items are random-access
  iterators. Normally, p points to the first element of
  a sequence and r to the first element beyond the sequence.

  Compiled and executed at DIKU as follows:
  % g++ -O3 selection.cc
  % a.out

*/

#include <algorithm> // std::min, std::nth_element, std::swap
#include <cassert> // macro assert
#include <cstdlib> // random, srandom, std::size_t, std::ptrdiff_t
#include <ctime> // std::clock, std::clock_t, macro CLOCKS_PER_SEC
#include <functional> // std::less
#include <iostream> // std::cout, std::endl
#include <iterator> // std::iterator_traits
#include <utility> // std::pair

template <class item>
void copy(item p, item q, item r) {
  item i;
  item j;
  for (i = p, j = r ; i < q; ++i, ++j)
    *j = *i;
}

template <class item, class ordering>
void insertion_sort(item p, item r, ordering less) {
  for (item j = p + 1; j < r; ++j) {
    typename std::iterator_traits<item>::value_type key = *j;
    item i = j - 1;
    while (i >= p && less(key, *i)) {
      *(i + 1) = *i;
      --i;
    }
    *(i + 1) = key;
  }
}

template <class item>
inline void insertion_sort(item p, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  insertion_sort(p, r, std::less<T>());
}

template <class item, class ordering>
void merge(item p, item q, item r, ordering less) {
  item s = new typename std::iterator_traits<item>::value_type[q - p];
  item t = s + (q - p);
  copy(p, q, s);
  item i;
  item j;
  item k;
  for (i = s, j = q, k = p; k < r; ) {
    if (less(*i,*j)) {
      *k = *i;
      ++i; ++k;
      if (i == t)
       break;
    }
    else {
      *k = *j;
      ++j; ++k;
      if (j == r) {
	copy(i, t, k);
	break;
      }
    }
  }
  delete[] s;
}

template <class item, class ordering>
void merge_sort(item p, item r, ordering less) {
  if (r - p > 1) {
    item q = p + (r - p) / 2;
    merge_sort(p, q, less);
    merge_sort(q, r, less);
    merge(p, q, r, less);
  }
}

template <class item>
inline void merge_sort(item p, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  merge_sort(p, r, std::less<T>());
}

template <class item, class ordering>
bool sorted(item p, item r, ordering less) {
  for (item q = p + 1; q < r; ++q) {
    if (less(*q, *(q - 1)))
      return false;
  }
  return true;
}

template <class item>
inline bool sorted(item p, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  return sorted(p, r, std::less<T>());
}

template <class item, class ordering>
void select_with_merge_sort(item p, item q, item r, ordering less) {
  merge_sort(p, r, less);
}

template <class item>
inline void select_with_merge_sort(item p, item q, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  merge_sort(p, r, std::less<T>());
}

template <class item, class ordering>
void select_with_STL_sort(item p, item q, item r, ordering less) {
  sort(p, r, less);
}

template <class item>
inline void select_with_STL_sort(item p, item q, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  sort(p, r, std::less<T>());
}

template <class item>
item pivot(item p, item r) {
  return p + (random() % (r - p));
}

template <class item, class ordering>
item partition(item p, item q, item r, ordering less) {
  typename std::iterator_traits<item>::value_type x = *q;
  std::swap(*p, *q);
  item i = p;
  item j = r;
  for (;;) {
    do i++; while (i < r && less(*i, x));
    do j--; while (less(x, *j));
    if (i < j) {
      std::swap(*i, *j);
    }
    else {
      std::swap(*p, *j);
      return j;
    }
  }
}

template <class item>
inline item partition(item p, item q, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  return partition(p, q, r, std::less<T>());
}

template <class item, class ordering>
bool partitioned(item p, item q, item r, ordering less) {
  for (item i = p; i < q; ++i) {
    if (less(*q, *i))
      return false;
  }
  for (item i = q + 1; i < r; ++i) {
    if (less(*i, *q))
      return false;
  }
  return true;
}

template <class item>
inline bool partitioned(item p, item q, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  return partitioned(p, q, r, std::less<T>());
}

template <class item, class ordering>
void randomized_select(item p, item q, item r, ordering less) {
  if (r - p > 1) {
    item p_prime = pivot(p, r);
    item q_prime = partition(p, p_prime, r, less);
    if (q < q_prime)
      randomized_select(p, q, q_prime, less);
    else if (q > q_prime)
      randomized_select(q_prime + 1, q, r, less);
  }
}

template <class item>
inline void randomized_select(item p, item q, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  return randomized_select(p, q, r, std::less<T>());
}

template <class item, class ordering>
std::pair<item, item> partition_3_way(item p, item q, item r, ordering less) {
  typename std::iterator_traits<item>::value_type v = *q;
  std::swap(*p, *q);
  item pa = p + 1;
  item pb = p + 1;
  item pc = r - 1;
  item pd = r - 1;
  for (;;) {
    while (pb <= pc && !less(v, *pb)) {
      if (less(*pb, v))
	pb += 1;
      else {
	std::swap(*pa, *pb);
	pa += 1;
	pb += 1;
      }
    }
    while (pb <= pc && !less(*pc, v)) {
      if (less(v, *pc))
	pc -= 1;
      else {
	std::swap(*pc, *pd);
	pd -= 1;
	pc -= 1;
      }
    }
    if (pb > pc)
      break;
    std::swap(*pb, *pc);
    pb += 1;
    pc -= 1;
  }
  std::ptrdiff_t smaller = pb - pa;
  std::ptrdiff_t i = std::min(pa - p, smaller);
  item k;
  item l;
  for (k = p, l = pb - i; l < pb; ++k, ++l)
    std::swap(*k, *l);
  std::ptrdiff_t larger = pd - pc;
  std::ptrdiff_t j = std::min(larger, r - pd - 1);
  for (k = pb, l = r - j; l < r; ++k, ++l)
    std::swap(*k, *l);
  return std::pair<item, item>(p + smaller, r - larger);
}

template <class item>
inline std::pair<item, item> partition_3_way(item p, item q, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  return partition_3_way(p, q, r, std::less<T>());
}

template <class item, class ordering, std::size_t gs>
void select(item p, item q, item r, ordering less) {
  if (r - p < 90)
    merge_sort(p, r, less);
  else {
    item i = p;
    item j;
    for (j = p; j < r - gs; ++i, j += gs) {
      insertion_sort(j, j + gs, less);
      std::swap(*i, *(j + gs / 2));
    }
    item p_prime = p + (i - p) / 2;
    select<item, ordering, gs>(p, p_prime, i, less);
    std::pair<item, item> z = partition_3_way(p, p_prime, r, less);
    if (q < z.first)
      select<item, ordering, gs>(p, q, z.first, less);
    else if (q >= z.second)
      select<item, ordering, gs>(z.second, q, r, less);
  }
}

template <class item, std::size_t gs>
inline void select(item p, item q, item r) {
  typedef typename std::iterator_traits<item>::value_type T;
  select<item, std::less<T>, gs>(p, q, r, std::less<T>());
}

void generate_int_sequence(int* p, int* r) {
   for (int* q = p; q < r; ++q)
	   *q = (int)q;
}

int main(int argc, char* argv[]) {
  const unsigned int repetitions = 2;
  const unsigned int n = atoi(argv[1]);
  int* a = new int[n];

  /*std::cerr << "Testing various selection algorithms..." << std::endl
            << "Select the k'th largest of the n random integers" << std::endl
            << "n: " << n << std::endl
            << "k: " << (n + 1)/2 << std::endl
            << "repetitions: " << repetitions << std::endl << std::endl;*/

  std::clock_t total_time;

  std::cout << "mergesort:";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    select_with_merge_sort(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << ""
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "STL sort:";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    select_with_STL_sort(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << ""
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "RND:";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    randomized_select(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << " "
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "Det (g = 03):";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    select<int*, 3>(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << " "
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "Det (g = 07):";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    select<int*, 7>(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << ""
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "Det (g = 11):";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    select<int*, 11>(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << ""
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "Det (g = 15):";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    select<int*, 15>(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << ""
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "Det (g = 19):";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    select<int*, 19>(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << ""
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "Det (g = 23):";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    select<int*, 23>(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << ""
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  std::cout << "STL sel:";
  srandom(4033);
  total_time = 0;
  for (unsigned int r = 0; r < repetitions; ++r) {
    generate_int_sequence(a, a + n);
    std::clock_t one_time = std::clock();
    std::nth_element(a, a + n / 2, a + n);
    one_time = std::clock() - one_time;
    total_time += one_time;
    assert(partitioned(a, a + n / 2, a + n));
  }
  std::cout << ""
            << float(total_time) / (float(repetitions) * float(CLOCKS_PER_SEC))
	    << std::endl;

  delete[] a;
}






