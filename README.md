# EIT
Algoritma rekonstruksi EIT berbasis Finite Element Model menggunakan library Eidors pada Octave/Matlab  
Website eidors: http://eidors3d.sourceforge.net/

Cara menggunakan EIDORS:
1. Unzip software eidors pada directory yang sudah ditentukan
2. Start Matlab dan run command run /path/to/eidors/startup.m
3. Coba jalankan program eidors_demo atau eidors_demo2 di folder EIT-master

# Simulasi
Folder ini berisi program simulasi rekonstruksi menggunakan metode Newton-Raphson (Regularisasi, Pseudoinverse, Singular Value Decomposition) dan Iterasi Landweber.  
Untuk mencoba simulasi, copy kan file-file simulasi kecuali ChooseCircle.m ke directory EIT-master\demo_complete.  
File ChooseCircle.m ini telah diupdate untuk dapat digunakan di octave yang tidak memiliki fitur untuk menentukan posisi dan ukuran lingkaran.  

# Program Rekonstruksi
Folder ini berisi program rekonstruksi EIT menggunakan data real tegangan bidang batas pada file excel.  
Untuk mencoba program rekonstruksi EIT ini, copy kan file-file program dan data excel ke directory EIT-master\demo_complete.

# Mesh Generator
Mesh generator menggunakan Netgen 6.1 yang memiliki requirement instalasi:
1. Python 3.7
2. Microsoft visual studio 2017 redistributable
3. EIDORS telah di-running

Lebih baik untuk melakukan instalasi Netgen pada directory EIDORS
Dokumentasi fungsi-fungsi dalam Netgen http://eidors3d.sourceforge.net/doc/index.html?eidors/meshing/netgen/call_netgen.html
Contoh pembuatan mesh terdapat pada folder Netgen

# Dokumentasi
Berisi dokumentasi progress penelitian EIT.  
