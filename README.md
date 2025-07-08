# slice-viewer
A set of tools for visualizing slices of 3D medical imaging data including MRI and CT.

## Notes
* This viewer depends on the NIfTI structure returned by the `nii_tool` function from `dicm2nii` repository. Many thanks to GitHub user `xiangruili` for the use of their repository (under the associated open licence).

## Initialization
* Clone this repository
```bash
git clone https://github.com/llawrence-main/slice-viewer.git
```
* Initialize the submodules:
```bash
git submodule init
git submodule update
```

## Usage
* Open MATLAB and run `example.m` to see an example of how to use the main function `view_slice` to visualize a slice of a 3D image.
* Check the comments at the top of each function for more information.
