#include "matmul.h"
#include <vector>
#include <thread>

void matmul_ref(const int* const matrixA, const int* const matrixB,
                int* const matrixC, const int n) {
  // You can assume matrixC is initialized with zero
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < n; k++)
        matrixC[i * n + j] += matrixA[i * n + k] * matrixB[k * n + j];
}

void matmul_optimized(const int* const matrixA, const int* const matrixB,
                      int* const matrixC, const int n) {
  // TODO: Implement your code
  
  std::vector<std::thread> threads;

  // REFERENCE: https://sgc109.tistory.com/130
  if (n == (n&-n) && n > 1024) {
    // when n is the power of 2
    const int n_threads = 128; // the matrix is big enough
    const int l_thread = n / n_threads; 
    // the number of columns of matrixB
    // that each thread needs to handle	

    for (int t = 0; t < n_threads; t++) {
      threads.push_back(std::thread([matrixA, matrixB, matrixC, n, l_thread, t, n_threads]{
	int* B_T = new int[l_thread * n];
	
	// Transpose a part of the matrix B
	for (int i = 0; i < n; i++) {
	  for (int j = 0; j < l_thread; j++) {
	    B_T[j * n + i] = matrixB[i * n + (t * l_thread + j)];
	  }
	}
        int e;
	for (int i = 0; i < n; i++) {
	  for (int j = 0; j < l_thread; j++) {
	    e = 0;
	    for (int k = 0; k < n; k++) {
	      e += matrixA[i * n + k] * B_T[j * n + k];
	    }
            matrixC[(i) * n + (j + (t*l_thread))] = e; 
	  }
	}
      }));
    }
  } else {
    int n_threads = 64;
    int l_thread = n / n_threads; // for convenience, the last thread is longer...
    // need to be optimzied further
    
    for (int t = 0; t < n_threads; t++) {
      threads.push_back(std::thread([matrixA, matrixB, matrixC, n, l_thread, t, n_threads]{
        const int c_optimal = 1024;
        int l_chunk;
        int n_chunks;

        if (n < c_optimal) {
          l_chunk = n;
	  n_chunks = 1;
        } else {
          l_chunk = c_optimal;
          n_chunks = n % l_chunk == 0 ? n / l_chunk : n / l_chunk + 1;
        }

        const int j_start = t * l_thread;
        const int j_end = t == n_threads - 1 ? n : j_start + l_thread;
        const int j_size = j_end - j_start;

        int* B_T = new int[j_size * n];

        for (int i = 0; i < n; i++) {
          for (int j = 0; j < j_size; j++) {
	    B_T[j * n + i] = matrixB[i * n + j + j_start];
	  }
        }

        int e;
        int k_start;
        int k_end;
        for (int c = 0; c < n_chunks; c++) {
          for (int i = 0; i < n; i++) {
      	    for (int j = 0; j < j_size; j++) {
	      e = 0;
	      k_start = c * l_chunk;
	      k_end = c == n_chunks - 1 ? n : k_start + l_chunk;
	      // the last thread will take a longer job if n is not divisible by l_chunk
	      for (int k = k_start; k < k_end; k++) {
	        e += matrixA[i * n + k] * B_T[j * n + k];
	      }
	      if (c == 0) {
	        matrixC[i * n + j + j_start] = e;
	      } else {
	        matrixC[i * n + j + j_start] += e;
	      }
	    }
	  }
        }
      }));
    }
  }

  for (auto& thread: threads) {
    thread.join();
  }
}
