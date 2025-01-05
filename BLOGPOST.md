# Towards the Standardization of Lattice-based Proof Systems through GPU Acceleration
#### *EMILE HREICH, 10.01.2025*

---

### Table of Contents
1. [Introduction](#introduction)
2. [Related Work](#related)
3. [Hardware Acceleration](#hardware)
4. [Design: Integrating ICICLE with lattirust-arithmetic](#lattirust)
5. [Implementation](#lova)
6. [Evaluation & Benchmarking](#applications)
7. [Concluding Remarks & Future Work](#conclusion)
8. [Acknowledgements](#acknowledgements)
9. [References](#references)

### Glossary
- **MPC**   : Multi-Party Computation
- **GPU**   : Graphics Processing Unit
- **FHE**   : Fully Homomorphic Encryption
- **TFHE**  : An open-source library for fully homomorphic encryption. [2](#FFHE) 
- **NTT**   : Number Theoretic Transform
- **SIS**   : Short Integer Solution
- **GPC**   : General-Purpose Computing
- **ASIC**  : Application-Specific Integrated Circuit
- **FPGA**  : Field-Programmable Gate Array
- **SIMD**  : Single Instruction, Multiple Data
- **AVX**   : Advanced Vector Extensions
- **SSE**   : Streaming SIMD Extensions
- **SIS**   : Short Integer Solution
- **MSIS**  : Module Short Integer Solution
- **CSR**   : Compressed Sparse Row
- **CSC**   : Compressed Sparse Column


---

### Introduction <a name="introduction"></a>
Lattices provide a versatile foundation for a broad spectrum of quantum-secure cryptographic schemes. The additional structural properties inherent in lattice-based assumptions offer a significant advantage over other post-quantum assumptions, such as those based on collision-resistant hash functions. This added structure allows for the design of cryptographic protocols that exhibit greater efficiency in their implementation. The operations used in lattice-based cryptography (and especially in succinct and zero-knowledge proofs) are essentially linear algebra operations over a ring (or similarly parallelizable computations).  Given the critical importance of performance for the adoption and standardization of such cryptosystems, this project explores the use of GPUs for advanced lattice-based cryptography, in particular by implementing a GPU backend for the arithmetics in the Lattirust library.

Lattirust is a Rust library designed to provide a secure and efficient foundation for lattice-based cryptography, with a focus on zero-knowledge proofs. The goal of Lattirust is to serve a role similar to arkworks for lattice-based cryptographic frameworks and lattigo for zero-knowledge and succinct proof systems. At its core, Lattirust provides a package, `lattirust-arithmetic`, for arithmetic operations used across the library. It implements power-of-two cyclotomic rings $\mathbb{Z}_q[X]/(X^{2^k}+1)$, number-theoretic transforms, and operations on matrices and vectors as well as various norms and challenge sets. The package also includes wrappers and helper traits for the nimue library, streamlining secure instantiations of the Fiat-Shamir transformation for lattice protocols. In terms of protocols, Lattirust includes an implementation of the LaBRADOR protocol, which delivers small proof sizes with linear verifier runtime, making it suitable for scenarios prioritizing compactness over verification speed. It also features lova, a lattice folding scheme designed for efficient recursive verification, enabling multiple computational relations to be folded into a single instance while maintaining controlled norm growth and reducing computational overhead.

To fully harness the computational capabilities of GPUs, this work focuses on accelerating key operations within the `lattirust-arithmetic` package, including matrix-matrix and matrix-vector operations, decomposition and recomposition processes, and number-theoretic transforms (NTT) for fast polynomial multiplication, and modular arithmetic. This post outlines a strategy for accelerating lattice-based cryptographic protocols on Nvidia GPUs by utilizing Ingonyama’s CUDA Backend. The primary goals are twofold: (1) achieving end-to-end integration with the ICICLE backend to enable the first deployment of the accelerated library and to evaluate the scope and implications of GPU acceleration, and (2) establishing both the theoretical and practical foundations for the eventual implementation of custom kernels tailored for further acceleration. This approach lays the groundwork for robust and scalable GPU-accelerated cryptographic solutions. At the time of writing, our contributions to Lattirust include:

- Conducting research into GPU acceleration for arithmetic operations implemented in Lattirust and providing a detailed specification for integrating the ICICLE GPU backend.  
- Developing the accelerated version of the `lattirust-arithmetic` package with GPU-accelerated support for critical operations such as NTT, decomposition, recomposition, and matrix-vector computations. This implementation is seamlessly integrated with Ingonyama's ICICLE backend.  
- Performing an experimental evaluation of the accelerated library implementation and presenting preliminary benchmark results to validate its performance improvements.  

### Related Work <a name="related"></a>

Research on accelerating cryptographic schemes using GPU architectures has gained substantial momentum. This is driven by the need to meet performance and scalability standards essential for the adoption and standardization of these schemes. The computational demands of modern security applications in privacy-preserving computation and post-quantum cryptography, make accelerating the operations that underpin these schemes critical to achieving practical performance and usability.

The available research on GPU acceleration in cryptography has primarily focused on Homomorphic Encryption and Secure Multi-Party Computation. These schemes are computationally intensive and require significant optimization to achieve practical performance. The use of GPUs in accelerating these schemes has shown promising results, with significant speedups achieved in various operations, such as basic Liner Algebra operations, Number Theoretic Transforms, and Modular Arithmetic. Several works on accelerating HE schemes exist in the literature. Early works in this field focus on evaluating the feasibility of GPU acceleration for basic computations on encrypted data. For instance, the work by [3](#CKKS) presents a GPU-accelerated implementation of the CKKS encryption scheme, a homomorphic encryption scheme that supports approximate arithmetic operations, for logistic regression over encrypted data. The authors demonstrate that while encrypted training remains slower than plaintext computation, their approach achieves comparable model accuracy, highlighting the potential for privacy-preserving ML applications. Their results suggest that FHE can support more complex ML applications, with performance improvements paving the way for broader privacy-preserving analytics on sensitive data.

More recent works focus on the performance comparison of different homomorphic encryption schemes accelerated on GPUs. The Phantom library [4](#Pantom) represents an advancement in leveraging GPUs to optimize word-wise HE schemes, specifically BGV, BFV, and CKKS, which support batch processing for privacy-preserving applications. This library introduces several innovations in GPU processing, such as kernel fusion, memory pooling, and efficient use of NTTs for polynomial multiplications. Phantom's benchmarks reveal substantial speed improvements over existing implementations, which reinforces the notion that GPU acceleration presents a viable solution for high-throughput applications.

Another contribution is REDsec [5](#REDsec), a framework that optimizes the execution of encrypted neural networks on GPUs. REDsec introduces the (RED)cuFHE library, which enhances GPU support for TFHE-based operations by enabling multi-GPU acceleration, leveled operations, and efficient handling of ciphertext noise through lazy bootstrapping techniques. REDsec’s optimizations enable faster encrypted inference on complex neural network models, demonstrating that neural network homomorphic operations can achieve manageable latency.

Both Phantom and REDsec mark a shift from feasibility studies toward practical, high-performance implementations of HE schemes on GPUs. Their innovations bring GPU-acceleration closer to meeting industry standards emphasizing the role of hardware acceleration in advancing complex cryptosystems' applicability to real-world scenarios. These works collectively underscore the trend towards standardizing GPU-accelerated operations, nonetheless, the application of GPUs to lattice-based cryptographic schemes in the context of proof-systems remains an underexplored area. 

Accelerating lattice-based proof systems using GPUs presents unique challenges compared to the GPU acceleration of traditional cryptosystems. Traditional cryptosystems, often rely on arithmetic operations that map well to GPU architectures, enabling efficient parallelization. In contrast, lattice-based cryptography involves complex mathematical structures and operations, such as high-dimensional polynomial multiplications, which are less straightforward to adapt and parallelize effectively on GPUs. One significant challenge in accelerating lattice-based cryptography on GPUs is the efficient implementation of polynomial multiplication, a core operation in these systems. NTTs are employed to perform these multiplications efficiently. However, implementing NTT on GPUs requires careful optimization to manage memory access patterns and parallel execution, as highlighted in [ADD REF](). Section III of this document provides a high level overview of the algorithms used for NTT acceleration on GPUs.

These inherent complexities continue to impede efficient GPU acceleration. This project seeks to address these challenges by exploring GPU acceleration strategies tailored for lattice-based cryptographic schemes. It focuses on critical operations underpinning these systems and evaluates the tools and techniques available to achieve efficient implementation.

### Hardware Acceleration <a name="hardware"></a>

The traction gained by GPU acceleration in cryptographic schemes stems not only from their inherent design for parallel processing but also from their ability to strike a practical balance between general-purpose computing platforms and specialized platforms like FPGAs and ASICs. GPUs are a versatile and cost-effective choice, suitable for a broad range of cryptographic applications. They boast shorter development cycles compared to FPGAs and ASICs, enabling faster prototyping and deployment, making them accessible to a wider range of users. Additionally, GPUs integrate seamlessly into existing systems, leveraging general-purpose computing infrastructure without necessitating specialized hardware modifications. This combination of adaptability, performance, and cost-efficiency positions GPUs as an advantageous middle ground for accelerating cryptographic workloads, including lattice-based cryptographic schemes.

#### SIMD on CPU vs. GPU

CPUs incorporate SIMD (Single Instruction, Multiple Data) capabilities through various instruction set extensions, such as Intel's SSE (Streaming SIMD Extensions) and AVX (Advanced Vector Extensions). These extensions enable CPUs to perform parallel operations on multiple data points simultaneously, enhancing performance in tasks like multimedia processing and scientific computations.However, the SIMD capabilities of CPUs are limited compared to GPUs. CPUs are optimized for general-purpose computing, emphasizing low-latency execution of single-threaded tasks. Their SIMD units typically feature narrower vector widths.Additionally, CPUs have fewer SIMD registers and limited memory bandwidth compared to GPUs, which are architected to handle high-throughput parallel workloads with thousands of concurrent threads. This architectural difference means that while CPUs can leverage SIMD for performance gains, they are less suited for tasks that require extensive parallelism, where GPUs excel.

#### GPU vs. FPGA

FPGAs are known for their power efficiency and ability to be customized for specific algorithms, often achieving higher performance per watt for tasks like modular arithmetic and hashing, common in cryptographic workloads. GPUs, by contrast, provide significant raw computational power and excel at parallel processing, making them well-suited for applications involving large-scale data parallelism. While GPUs often outperform FPGAs in terms of absolute performance for highly parallel tasks, FPGAs are better suited for fine-grained parallelism. Developing for FPGAs requires expertise in hardware description languages and complex synthesis processes, whereas GPUs benefit from mature programming frameworks like CUDA and ICICLE, streamlining development and prototyping.

### Design: Integrating ICICLE with lattirust-arithmetic <a name="lattirust"></a>

The *lattirust-arithmetic* package is responsible for handling a variety of algebraic operations central to the library’s GPU-accelerated functionality. We focus on five of its core operations which are bottlenecks in the lattice-based cryptographic protocols that Lattirust supports:

1. **Decomposition**  
2. **Recomposition**  
3. **Number Theoretic Transforms (NTTs)**  
4. **Inner Products**  
5. **Sparse Matrix-Vector Multiplication**

To support these operations efficiently, several key parameters must be considered:

- **Degree of the cyclotomic polynomial ($\boldsymbol{D}$):** Common choices for $(D)$ are 64, 128, and 256.  
- **Field modulus ($\boldsymbol{q})$:** This is chosen as an NTT-friendly prime appropriate for the selected $(D)$. It typically ranges between 30 and 60 bits, and its special structure does not compromise security.  
- **Vector/matrix dimensions:** These may range from large witness lengths (on the order of $(2^{20})$ to $(2^{30})$ ) down to more modest sizes (128 to 1024).  
- **Decomposition basis $(\boldsymbol{b}$):** Commonly set to either $(b = 2)$ or $(b = \sqrt{q})$.


#### Decomposition & Recomposition

The decomposition operations is handled in the `decompose_balanced` function of the package. It is implemented for rings and field elements but its functionality also extends to vectors and matrices, where each element of a matrix can be decomposed according to the same principles. The function decomposes field elements with respect to a basis `b` into an equivalent representation of multiple field elements with smaller norm. We consider a balanced decomposition, where each element of the decomposition lies in the range $[-b/2, b/2]$. The basis `b` must be greater than 1. The current implementation requires `b` to be even for simplicity. 

The implementation iteratively extracts the remainder (`rem`) of the current value divided by the basis `b`. This remainder represents the next component in the decomposition sequence. If `rem` falls within the $([-b/2, b/2])$ range, it is appended to the decomposition, and the current value is reduced by dividing it by `b`. If `rem` falls outside this range, adjustments are made to bring it back within bounds by adding or subtracting `b`, followed by computing and applying a "carry" to adjust the current value. This step ensures each component of the decomposition stays within the balanced range. The loop terminates once the current value reaches zero, meaning the entire element has been decomposed. If a padding size is specified, the function pads the decomposition to this length with zeros. Otherwise, it removes any trailing zeros that don’t affect the original value.

The recomposition function complements decomposition, enabling the reassembly of previously decomposed elements back into their original form. The `recompose` function takes a vector of `FinitefieldWrapper` elements and reconstructs the original finite field element. This function essentially sums the values weighted by their respective powers of `b`, reversing the decomposition process.

This functionality also extends to matrices. Two additional functions are implemented to that end. The `recompose_matrix` function is used for matrices where each row has been decomposed into multiple columns using a basis. If $A$ is a $m\times n$-matrix decomposed into a $m\times nk$ matrix $\tilde{A}$, then $A = \tilde{A}G$ for the gadget matrix $G = g \otimes I$ where $g=\left[\begin{array}{c}1 \\b\\\vdots\\b^k\end{array}\right]$ and $I$ is an identity matrix of suitable dimensions. The `recompose_left_right_symmetric_matrix` function performs recomposition on symmetric matrices (in particular, matrices $M = A^\top A$), where the structure is $G^T * mat * G$, and `G` represents the gadget matrix. The function iterates through pairs of row and column indices, recomposing each element based on the `powers_of_basis`. For each symmetric matrix element, it multiplies the decomposed component by the corresponding product of powers, summing these products to reassemble the original element. The final result is an `n x n` symmetric matrix.

Overall, these operations involve loop-carried dependencies and a relatively complex control flow inside the loop, making them a poor fit for GPU acceleration. GPUs are generally optimized for data-parallel workloads rather than speeding up a single, serial chain of computation, so the core decomposition and recomposition operations do not benefit much from a GPU approach. However, when extending these operations to vectors and matrices (via decompose_balanced_vec, decompose_matrix and their complementary recomposition functions), they are simply applied element-wise. In this context, minimal acceleration can be achieved on CPU, leveraging SIMD instructions for example, by employing batched decomposition and recomposition, where multiple elements from the vector or matrix are processed in parallel.

#### Number Theoretic Transforms

The Number Theoretic Transform is a generalization of the Discrete Fourier Transform over finite fields, where the field is defined by a prime number $q$ and a primitive root of unity $\omega$. The NTT is used to convert a polynomial from the time domain to the frequency domain, facilitating efficient polynomial multiplication through pointwise multiplication in the frequency domain. Efficient NTT implementations on GPUs reduce the time complexity of polynomial multiplications from quadratic $O(n^2)$ to quasi-linear $O(n\log{}n)$.

##### Forward NTT
The forward NTT transforms a polynomial $f(x)= f_0+f_1x+...+f_{n-1}x^{n-1}$ into its NTT representation, evaluating the polynomial at roots of unity. For a polynomial $f(x)$ of degree $n−1$, the forward NTT computes the following transformation:

$F(k) = \sum_{j=0}^{n-1}f_j\omega^{j{}k}{}\mod{}q$ for $k = 0, 1, \dots, n-1$

Where:
- $q$ is a large prime modulus,
- $\omega$ is an n-th primitive root of unity such that $\omega{}n≡1 \mod q$

##### Inverse NTT

Polynomial multiplication is performed via component-wise multiplication of the transformed coefficients in the NTT space. To retrieve the original polynomial from the NTT-transformed representation, the inverse NTT is used. It’s given by:

$f_j = n^{-1}\sum_{k=0}^{n-1}F(k)\omega^{-j{}k}{}\mod{}q$ for $j = 0, 1, \dots, n-1$

Where $n^{-1}$ is the modular inverse of $n$ modulo $q$. This inverse transform restores the original polynomial coefficients after pointwise multiplication in the NTT domain.

##### Accelerating NTT on GPUs

There exists two specialized algorithms for the acceleration of NTT: the Cooley-Tukey Butterfly Algorithm [ADD REF]() and the Four-Step NTT Algorithm [ADD REF](). Each algorithm is tailored for different use cases based on polynomial degree, available hardware, and the nature of the cryptographic task.

###### Cooley-Tukey Butterfly Algorithm
This algorithm represents the classical approach for computing the NTT. It is particularly well-suited for smaller polynomials or scenarios where simplicity and low overhead are essential. The process involves iterative, stage-wise computations of polynomial coefficients using modular arithmetic. It is well-suited for problems where the degree n of the polynomial is relatively small (e.g.,n ≤ $2^{16}$) and the memory or thread count available on the GPU is sufficient to handle all the coefficients in parallel.

The first step is the **bit-reversal permutation**, where the input vector is rearranged into bit-reversed order. This preprocessing ensures that subsequent computations can be performed efficiently in place, simplifying the overall algorithm.

The second step involves **butterfly operations**, which constitute the core of the computation. At each stage, coefficients are combined using modular addition and multiplication by powers of a primitive root of unity $\omega$. The results are then stored back in the input vector, progressively halving the distance between paired indices at each stage. This systematic approach allows for efficient and low-overhead computation of the NTT.

###### Four-Step NTT Algorithm

The Four-Step NTT Algorithm is designed to address the memory access challenges associated with larger polynomials, particularly when $ n \geq 2^{20} $. This algorithm employs a divide-and-conquer strategy to split the problem into smaller, more manageable subproblems, thereby reducing memory access latency and improving cache locality.

The algorithm begins with a matrix transposition, where the input polynomial is reshaped into a 2D matrix with dimensions $n_1 \times n_2 $. This restructuring improves memory access patterns by enabling smaller, independent NTTs to be performed on each row and column of the matrix. Following this, independent $ n_1$-point NTTs are applied to each row of the matrix. This step is fully parallelizable, allowing each row to be processed by a separate set of GPU threads, maximizing computational throughput.

After the row-wise NTTs, each element of the matrix undergoes a twiddle factor multiplication. This step aligns the results correctly by adjusting for the separation into smaller transforms. The algorithm then applies a column NTT, where $n_2$-point NTTs are performed on each column of the matrix. This completes the transformation of the polynomial coefficients. Finally, a second matrix transposition restores the original layout of the polynomial coefficients, ensuring the output matches the expected format.

###### ICICLE CUDA Backend

Ingonyama's ICICLE CUDA Backend provides an implementation that supports two algorithms for NTT acceleration: `radix-2` and `mixed-radix`. These algorithms are variations of the Cooley-Tukey Butterfly Algorithm, optimized for finite fields. The choice between them depends on the degree of the polynomial and the available hardware resources, as each has distinct computational and memory requirements.

The Radix-2 NTT algorithm is designed for input sequences whose lengths are powers of two, leveraging a divide-and-conquer approach to simplify computations. The input sequence is recursively divided into smaller subsequences of even-indexed and odd-indexed elements. This process halves the problem size at each stage. The core of the computation lies in the butterfly operation, where pairs of elements are combined using a modular arithmetic formula. Precomputed twiddle factors, which are roots of unity in the finite field, are used to optimize these operations. A final step reorders the output from a bit-reversed sequence into its natural order.

The butterfly operation for Radix-2 is:


$X_k = (A_k + B_k \cdot W_k) \mod p$

where:
- $( A_k )$: Element from the even-indexed subset.
- $( B_k )$: Element from the odd-indexed subset.
- $( W_k )$: Twiddle factor (precomputed root of unity).
- $( p )$: Prime modulus.

In contrast, the Mixed-Radix NTT algorithm extends Radix-2 to accommodate input sequences of arbitrary composite lengths. This flexibility is particularly advantageous for polynomials with degrees that are not powers of two. The Mixed-Radix approach decomposes input sequences into blocks of sizes such as 16, 32, or 64, depending on the factorization of the sequence length. Unlike Radix-2, which strictly divides into size-2 subsequences, Mixed-Radix utilizes a variety of radices, optimizing the number of stages required for larger inputs. The decomposition strategy allows for a more efficient use of memory and computational resources.

The generalized butterfly operation for Mixed-Radix is:

$X_{k, r} = \sum_{j=0}^{r-1} (A_{j,k} \cdot W_{j,k}) \mod p$

where:
- $( X_{k, r})$: Output for the $k$-th set of inputs.
- $( A_{j,k} )$: $j$-th input element for the $k$-th operation.
- $( W_{j,k} )$: Twiddle factor specific to the radix.
- $( p )$: Prime modulus.

After performing all butterfly operations, Mixed-Radix recombines the results into a single output sequence. Reordering in Mixed-Radix is more complex than in Radix-2, often requiring digit-reversal permutations due to the varied sizes of the sub-transforms.

The Radix-2 algorithm is best suited for input sizes that are powers of two, offering simplicity and computational efficiency. Meanwhile, the Mixed-Radix algorithm is ideal for inputs with composite lengths, providing greater flexibility and optimized performance for larger polynomials. 

#### Inner Products & Sparse Matrix-Vector Multiplication

The inner product operation computes the sum of element-wise products between two vectors in the field. The implementation computes the inner products between pairs of vectors, producing a scalar. `inner_products_serial` computes the inner product of each pair of vectors in a given vector set, returning the results as a symmetric matrix, while `inner_products` does the same thing with some basic parallelism. The `inner_products_mat` function takes a matrix $A$ and computes the symmetric matrix $A^\top A$ (which we store as an upper triangular matrix).  

The `sparse_matrix_vec_mult` function multiplies a sparse matrix (in Compressed Sparse Column (CSC) format) by a dense vector.

##### Accelerating Matrix Operations on GPUs

For general efficient matrix multiplication on GPU, we suggest to follow the Nvidia cuTLASS tiling structure. It decomposes the computation into a hierarchy of thread block tiles, warp tiles  and thread tiles, effectively mirroring the CUDA programming model. The cuTLASS library has supported operations on sparse matrices since the 2.3 release, making it a viable option for managing the sparse-dense operations. Alternatively, the cuSPARSE library can be utilized, which, like cuBLAS (closed source), provides a specialized set of kernels designed explicitly for handling operations on sparse matrices. 

###### ICICLE CUDA Backend

As of now, ICICLE supports only matrix transposition and lacks built-in matrix multiplication capabilities. However, in the lattirust framework, matrix-vector multiplication can still be optimized using ICICLE's existing vector multiplication functions. This is achieved by decomposing the matrix-vector multiplication into multiple vector-vector multiplications, where each row of the matrix is multiplied by the vector.

While this method is straightforward, it doesn't account for the matrix's sparsity or its in-memory storage format, which can lead to inefficiencies. Sparse matrices, characterized by a majority of zero elements, benefit from specialized storage formats that enhance computational performance and reduce memory usage.

**Compressed Sparse Row (CSR) Format:**

In CSR format, a matrix is represented using three arrays:

1. **Values:** Contains all the non-zero elements of the matrix.

2. **Column Indices:** Stores the column index corresponding to each value.

3. **Row Pointers:** An array where each entry indicates the starting position of a row in the 'Values' and 'Column Indices' arrays.

This structure is particularly efficient for row-wise operations, such as iterating over elements in a specific row, making it advantageous for algorithms that process matrices row by row.

**Compressed Sparse Column (CSC) Format:**

Similarly, the CSC format uses three arrays:

1. **Values:** Contains all the non-zero elements.

2. **Row Indices:** Stores the row index for each value.

3. **Column Pointers:** An array where each entry points to the starting position of a column in the 'Values' and 'Row Indices' arrays.

CSC is more efficient for column-wise operations, such as accessing or modifying elements in a particular column, and is beneficial when algorithms process matrices column by column.


### Implementation <a name="lova"></a>

The implementation of GPU-accelerated operations in the lattirust-arithmetic package is specifically designed to leverage the ICICLE framework's Rust bindings for the CUDA backend, optimizing performance on Nvidia GPUs. Ingonyama is actively broadening hardware support for their framework, with plans to include a Metal/Apple Silicon backend. This expansion aims to provide an efficient solution for client-side cryptographic operations. Further details on this development will be discussed in the final section of this post.

Lattirust handles the core arithmetic operations within the lattirust-arithmetic package. The implementation is designed to ensure that the library's broader functionalities remain agnostic to the underlying hardware and acceleration mechanisms. Consequently, the only package directly involved is lattirust-arithmetic, and the sole added dependency is the ICICLE framework. This approach allows for a clean separation of concerns, allowing for seamless replacement of the accelerated operations with alternative implementations in the future, if necessary. Additionally, we ensure that the interfaces in the library's original implementation remain unchanged, preserving compatibility with existing protocol implementations. Furthermore, the operations are designed to gracefully fall back to the CPU implementation in the event of any issues during GPU execution.

To manage the GPU-accelerated operations, we introduce a new module, `gpu_context`. Its primary role is to handle the initialization, configuration, and data transfer between the host (CPU) and the device (GPU). Additionally, it provides utility functions for efficient conversion between field elements used in cryptographic computations. The module acts as a bridge between the cryptographic algorithms implemented in lattirust and the GPU-accelerated backend provided by ICICLE. By abstracting GPU-specific details and providing fallback mechanisms, it ensures that the library remains robust and hardware-agnostic. Its clean separation of responsibilities simplifies debugging, extensibility, and maintenance.

#### Initialization and Backend Configuration
- `load_backend`attempts to load the ICICLE CUDA backend from the specified environment variable which locates the backend on the device. 
- `try_load_and_set_GPU_backend_device` attempts to register the GPU device and sets it as the active computation target. If GPU initialization fails, it falls back to CPU implementation.
- `Error Management` is handled with an enum that provides a categorization of possible issues that may arise such as:
    - Missing GPU Backend
    - Device unavailability
    - Data transfer failures

#### Data Conversion Utilities
To this day, ICICLE only supports operations on algebraic fields, whereas many operations in Lattirust are implemented over Rings. Further, the only compatible field in ICICLE is Babybear. Discussions are underway to expand ICICLE's capabilities to include fields tailored to Lattirust's  protocols specific requirements.

For **Labrador**, two fields from the set {Babybear, KoalaBear, TeddyBear} are utilized. Labrador constructs a SNARK directly, benefiting from fields with approximately 64-bit sizes for efficient computations. In R1CS reduction, the size of the smallest prime factor of the modulus $q$ is crucial. Employing two Bear subfields is advantageous, potentially resulting in smaller proof sizes compared to the Labrador split.

**Greyhound** requires a single ~32-bit prime $p$ where $p \equiv 5 \mod 8$. However, NTTs cannot be performed directly in this field. The approach involves selecting 12-14 primes that support NTTs, similar to those used in Gregor's implementation, and employing the Bernstein-Sorenson technique to execute NTTs within this CRT basis. Integrating these components requires additional effort to ensure seamless functionality.

For **Lova**, there are no specific modulus constraints, and it utilizes $q = 2^{64}$. While Lova and other unstructured constructions are significant, they are currently less critical compared to the highly optimized structured constructions.

In this work, we provide a preliminary implementation that uses the Babybear field for GPU-accelerated operations. The gpu_context module includes functions to convert between the Babybear field and the Ring field used in Lattirust, ensuring the data is correctly formatted for GPU execution with the ICICLE CUDA backend. Two types of conversion functions are implemented. The first, convert_to_babybear and convert_to_zq, handles the conversion of field elements defined in Lattirust to Babybear elements and vice versa. The second, convert_poly_vector_to_babybear and convert_babybear_to_poly_vector, manages the conversion of PolyRing elements in Lattirust to Babybear elements and back.

##### `convert_to_babybear`

Converts an array of `Zq<Q>` elements into a `Vec<BabybearField>`.

**Implementation:**  
1. The function iterates over the input array of `Zq<Q>` elements.
2. Each `Zq<Q>` element is converted to its big integer representation using `.into_bigint()`. 
3. The first 32 bits of this big integer are extracted (`.0[0] as u32`) and used to construct a `BabybearField` element.
4. The result is collected into a `Vec<BabybearField>` and returned.

---

##### `convert_poly_vector_to_babybear`

Converts a polynomial vector from a `PolyRing` representation into a `Vec<BabybearField>`.

**Implementation:**  
1. The polynomial coefficients are flattened into a single vector using `P::flattened(poly_vec)`.
2. Each coefficient (a base-ring element) is converted into its unsigned representative (`UnsignedRepresentative`).
3. The unsigned representative is then converted into a `BabybearField` using its first 32 bits (`usr.0 as u32`).
4. The resulting `BabybearField` elements are collected into a `Vec<BabybearField>`.

---

##### `convert_babybear_to_poly_vector`

Converts a vector of `BabybearField` elements back into a `Vector<P>` in the `PolyRing` format.

**Implementation:**  
1. Each `BabybearField` element is converted back to its base-ring representation using `P::BaseRing::from_bytes()` and its byte representation (`to_bytes_le()`).
2. These base-ring elements are grouped into chunks corresponding to the dimension of the polynomials (`P::dimension()`).
3. Each chunk is used to reconstruct a polynomial using `P::from()`.
4. The reconstructed polynomials are collected into a `Vector<P>`.

---

##### `convert_to_Zq`

Converts a vector of `BabybearField` elements into an array of `Zq<Q>` elements.

**Implementation:**  
1. The function first validates that the input vector's length matches the expected output array length.
2. Each `BabybearField` element is converted into its byte representation using `.to_bytes_le()`.
3. The raw value is reconstructed from these bytes using a bitwise operation (`(byte as u64) << (index * 8)`).
4. A new `Zq<Q>` element is created using the reconstructed value.
5. The converted elements are collected into an array and returned as a `Result`.

#### Memory Management

The `gpu_context` module manages memory allocation and data transfers as required for GPU-accelerated operations. Many ICICLE functions take a configuration parameter that specifies the location of the data (host or device). This approach offers two significant advantages. First, it simplifies development for applicable operations by making the programming process largely oblivious to data transfers—though this abstraction is not without risks. Second, it allows for the optimization of data transfer, ensuring efficient resource utilization and reduced overhead. In fact, in most lattice-based protocols, a matrix $A$ is typically decomposed into $\tilde{A}$, and then immediately committed to as $C = KA$  (lattice commitments is a dense-dense matrix-matrix multiplication). It would thus be advantageous to keep this data on GPU between these operations. Similarly, at a later stage of the protocol, a matrix $\tilde{A}$ is recomposed as $\tilde{A}G$, and then the L2 norm of its column vectors is checked (using inner products). 

#### NTT Context Management and Acceleration

The context module also provides essential functions to initialize and manage the NTT context on the GPU, enabling efficient polynomial transformations.

The `init_ntt_context_on_device` function is responsible for allocating memory on the GPU and setting up the NTT domain using a specified primitive root. Domain initialization is performed once per GPU session because twiddle factors—roots of unity in the finite field—are cached to save runtime. The current implementation uses ICICLE's default NTT configuration, as this is a preliminary implementation. This configuration allows customization of parameters such as batch size, column batch computation, and the order of inputs and outputs. Future work will focus on benchmarking with custom settings tailored to adapted fields.

The `gpu_ntt_acceleration` function handles both forward and inverse NTT operations. The direction is specified by the `dir` parameter in ICICLE's function signature. The process is as follows:
1. The function attempts to load and set the GPU backend using `try_load_and_set_GPU_backend_device`. If the GPU backend is unavailable, it falls back to a CPU implementation.
2. It initializes the GPU's NTT context via `init_ntt_context_on_device`, ensuring that input and output memory are allocated and prepared for computation.
3. Input data is converted from `Zq<Q>` elements to `BabybearField` elements using `convert_to_babybear` and transferred from host memory to GPU memory using `copy_from_host`.
4. The NTT operation is performed on the GPU with ICICLE's default NTT configuration.
5. The results are copied back to host memory using `copy_to_host` and converted from `BabybearField` elements back to `Zq<Q>` elements using `convert_to_Zq`.

If GPU acceleration is unavailable or disabled at compile time, the `gpu_ntt_acceleration` function defaults to a CPU-based implementation. While this approach is less efficient, it ensures that the system remains functional and reliable.

#### Vector Operations Acceleration

<>

### Evaluation & Benchmarking <a name="applications"></a>

The GPU we use in this work is the Nvidia 3090 with the following specifications:
- GPU Codename : GA102
- GPU Architecture : Ampere
- SMs : 82
- Cuda Cores / SM	: 128
- Cuda Cores / GPU : 10496
- Tensor Cores / SM : 4
- Tensor Cores / GPU : 328
- GPU Boost Clock (MHz) : 1695
- CUDA Compute Capability	: 8.6
- Peak FP32 TFLOPS : 35.6
- Peak FP16 TFLOPS : 35.6
- Memory Interface : 384-bit
- Memory Clock (Data Rate) : 19.5 Gbps
- Memory Bandwidth : 936 GB/sec
- L1 Data Cache	: 10496 KB
- L2 Cache Size : 6144 KB
- Register File Size : 20992 KB

### Concluding Remarks & Future Work <a name="conclusion"></a>

This project marks an advancement in accelerating Lattirust, a library dedicated to safe and efficient lattice-based cryptography. We have pinpointed the operations most crucial for acceleration and devised a comprehensive strategy for their implementation on GPU. The integration of the library with the ICICLE framework has been particularly advantageous during the prototyping and development phases. This framework, which conveniently includes Rust bindings, has proven very useful in navigating and adapting the Lattirust ecosystem. It has facilitated the identification of optimization opportunities, informed software engineering decisions, and, most importantly, the construction of preliminary benchmarks for the accelerated operations. Consequently, it has significantly reduced the time required to develop a functional version of the accelerated library.

Our implementation of GPU-accelerated operations highlights the feasibility and significant performance gains achievable in lattice-based cryptographic protocols through GPU acceleration. However, ICICLE is not yet the definitive solution for accelerating a library like Lattirust. Its current development is still in its early stages, supporting a limited range of algebraic structures and operations, with a primary focus on zero-knowledge proofs. Additionally, while the CPU implementation of ICICLE is fully open source, the CUDA backend is currently designated as "delayed open-source." Ingonyama intends to release the backend under an MIT license once the integration reaches a stable state. This limitation restricts the ability to customize and fine-tune the GPU kernels for specific use cases, which could otherwise enhance performance further. Despite these constraints, we have worked closely with Ingonyama’s team to mitigate the impact of these limitations, ensuring meaningful progress in adapting ICICLE for our requirements.

An area where ICICLE proves particularly useful is in managing complex algebraic structures. In the context of GPUs, the infrastructure required for these structures can vary significantly based on the configurations. Implementing cryptographic algorithms across different rings and fields necessitates customized approaches to address the distinct arithmetic properties and operational requirements of each algebraic structure. For finite fields, operations typically involve modular arithmetic with a prime modulus. Since GPUs do not natively support modular arithmetic, efficient implementation strategies are necessary. When working with rings, such as integers modulo n or polynomial rings, the arithmetic operations often become more complex. They may involve more complex arithmetic, such as polynomial multiplication or operations with composite moduli. Implementing these operations efficiently on GPUs requires specialized algorithms that can exploit parallelism. Furthermore, accross all implementations, memory bandwidth & latency, as well as, precision and numerical stability are critical considerations to take into account. The idea of Ingonyama's ICICLE framework is to provide a unified interface for these operations, allowing developers to focus on the cryptographic protocols rather than the low-level GPU optimizations.Besides, as mentioned in previous sections, Ingonyama plans to include a Metal/Apple Silicon backend in the future [ADD REF](), which would be beneficial on many aspects. Unlike CUDA, which relies on discrete memory architectures, Metal's unified memory model allows both the CPU and GPU to access the same memory pool without requiring costly data transfers. This eliminates the significant overhead associated with repeated CPU-to-GPU and GPU-to-CPU transfers, a major bottleneck in CUDA workflows, particularly when CPU computation is needed as a fallback. The unified memory advantage in Metal drastically improves performance for mixed workloads by reducing latency and simplifying memory management. For example, in workflows where data must be frequently accessed or modified by both the CPU and GPU, Metal avoids the need for explicit memory copying, resulting in near-zero transfer costs. This efficiency was highlighted in performance benchmarks performed by Ingonyama [ADD REF](), where Metal consistently outperformed CUDA, especially in tasks involving iterative data exchanges. Additionally, Metal's unified memory model aligns well with the requirements of Izero-knowledge proof computations, which often involve heavy interaction between different computational components.

In the perspective of improving and maintaining a high-quality infrastructure for the acceleration of Lattirust, we foresee several areas for future exploration and improvement. One area is to continue the close collaboration with Ingonyama to ensure that they offer the right infrastructure and implementations for Lattirust's specific needs. This collaboration entails providing detailed feedback on the limitations and bottlenecks of the current system, as well as offering precise specifications for new features or optimizations. On one hand, engaging in such a dialogue would help Ingonyama refine their frameworks and ensure that they are aligned with the practical requirements of lattice-based cryptographic acceleration, on the other hand, this collaboration brings Lattirust one step closer to standardization. Another critical area for future work relates to the quality and reliability of GPU-accelerated software. Subtle synchronization bugs in GPU programs can lead to intermittent failures, jeopardizing the reliability of cryptographic computations. Additionally, sub-optimal programming practices can result in performance bugs, causing underutilization of GPUs, which are inherently resource-intensive and power-hungry. Future work will involve a thorough analysis of synchronization mechanisms and a systematic approach to identify and mitigate performance bottlenecks in GPU programs. This could include the exploration of a tool named iGUARD, which can pinpoint synchronization bugs to programmers [ADD REF]() and ScopeAdvice, which can help programmers find performance bugs due to over and redundant synchronization in GPU programs [ADD REF]().

### Acknowledgements <a name="acknowledgements"></a>

Many thanks to Christian Knabenhans for the great guidance and discussions throughout the project.

### References <a name="references"></a>

```
@article{10428046,
      author={Yang, Hao and Shen, Shiyu and Dai, Wangchen and Zhou, Lu and Liu, Zhe and Zhao, Yunlei},
      journal={IEEE Transactions on Dependable and Secure Computing}, 
      title={Phantom: A CUDA-Accelerated Word-Wise Homomorphic Encryption Library}, 
      year={2024},
      volume={},
      number={},
      pages={1-12},
      doi={10.1109/TDSC.2024.3363900}
}
```

<a name="FFHE"></a>
```
@misc{cryptoeprint:2016/870,
      author = {Ilaria Chillotti and Nicolas Gama and Mariya Georgieva and Malika Izabachène},
      title = {Faster Fully Homomorphic Encryption: Bootstrapping in less than 0.1 Seconds},
      howpublished = {Cryptology {ePrint} Archive, Paper 2016/870},
      year = {2016},
      url = {https://eprint.iacr.org/2016/870}
}
```

<a name="CKKS"></a>
```
@misc{cryptoeprint:2018/462,
      author = {Hao Chen and Ran Gilad-Bachrach and Kyoohyung Han and Zhicong Huang and Amir Jalali and Kim Laine and Kristin Lauter},
      title = {Logistic regression over encrypted data from fully homomorphic encryption},
      howpublished = {Cryptology {ePrint} Archive, Paper 2018/462},
      year = {2018},
      url = {https://eprint.iacr.org/2018/462}
}
```

<a name="Phantom"></a>
```
@ARTICLE{10428046,
      author={Yang, Hao and Shen, Shiyu and Dai, Wangchen and Zhou, Lu and Liu, Zhe and Zhao, Yunlei},
      journal={IEEE Transactions on Dependable and Secure Computing}, 
      title={Phantom: A CUDA-Accelerated Word-Wise Homomorphic Encryption Library}, 
      year={2024},
      volume={21},
      number={5},
      pages={4895-4906},
      keywords={Graphics processing units;Homomorphic encryption;Finite element analysis;Phantoms;Libraries;Benchmark testing;Arithmetic;BFV;BGV;CKKS;GPU acceleration;homomorphic encryption},
      doi={10.1109/TDSC.2024.3363900}}

```

<a name="REDsec"></a>
```
@misc{cryptoeprint:2021/1100,
      author = {Lars Folkerts and Charles Gouert and Nektarios Georgios Tsoutsos},
      title = {{REDsec}: Running Encrypted Discretized Neural Networks in Seconds},
      howpublished = {Cryptology {ePrint} Archive, Paper 2021/1100},
      year = {2021},
      url = {https://eprint.iacr.org/2021/1100}
}
```


