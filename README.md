# CG-DNA-model
# A coarse-grained DNA model to study protein-DNA interactions and liquid-liquid phase separation: Application for mechanism and regulation beyond the Nucleosome
**Authors:** Utkarsh Kapoor [1], Young C. Kim [2], Jeetain Mittal [1] <br>
[1] Artie McFerrin Department of Chemical Engineering, Texas A&M University, College Station, Texas 78743, United States <br>
[2] Center for Materials Physics and Technology, Naval Research Laboratory, Washington, District of Columbia


This directory contains the source code required to run simulations of the 2-bead and 3-bead coarse-grained DNA models within the LAMMPS MD package. The directory also contains sample input files for running simulation of a DNA duplex at a given temperature. These coarse-grained models are highly-robust and are capable of studying processes such as hybridization and DNA protein interactions. For details about the model itself and it's development, please see Kapoor et al. (DOI: )

### 1. Source Codes
CG simulations were performed using the LAMMPS molecular dynamics simulations package (Oct 2020 version), in which CG-DNA-model codes have been implemented. To run these codes successfully one would require LAMMPS Oct 2020 package installed with the files within 1.LAMMPS_subroutines folder.

### 2. Examples
This directory contains sample input files for running simulation of a DNA duplex (S1: 5´-GCGTCATACAGTGC-3´ and its complement S2: 5´-GCACTGTATGACGC-3´) using 2-bead and 3-bead CG-DNA models.

#### LAMMPS Versions
The CG-DNA-model has not been tested with versions of LAMMPS other than Oct 2020. We welcome pull requests that update CG-DNA-model to be compatable with more recent versions of LAMMPS. In the meantime however, it is recommended that you use the Oct 2020 version of LAMMPS for your simulations. If you have any questions and concerns associated with bugs (not with general CG-DNA-model compilation or LAMMPS issues) please feel free to open a GitHub issue.

#### LICENSE AND DISCLAIMER
Redistribution and use of this software in source and binary forms, with or without modification, are permitted provided that this statement and the following disclaimer are retained. THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
