# Phantom Asphere

This is a MATLAB code to realize the Phantom Asphere method based on CODEV. The Phantom Asphere is a lens design technique where an asphere is added to the system to guide the designer in determining where to add a new element. The workflow is as follows:

1. **Add an Asphere**: Add an asphere to a singlet in the system. The designer can experiment with arbitrary singlets in the system to determine which one improves the performance the most.

2. **Replace with a Singlet**: Based on the asphere coefficients, replace the aspheric surface with a singlet (with a certain shape determined by the asphere coefficients).

3. **Optimize the System**: Re-optimize the system.

This method is very useful as it not only tells the designer the location of the new element but also its shape. However, determining the shape visually can be challenging. 

This script uses MATLAB's Global Optimization Toolbox to search for the surface shape that provides a similar aberration contribution as the asphere.

---

## How to Use

The purpose of this code is for step 2 mentioned above, directly tell the lens shape from the asphere coefficient is hard so I wrote a MATLAB code to address this problem:

1. **Open the Script**: Open `phantom_asphere_ver3.m`.

2. **Set the Surface Number**: Change the surface number in the script to the surface where you want to implement the Phantom Asphere.

3. **Run the Code**: Execute the script. The code will output six parameters in the following order:


- **R11**, **R12**: Radii of the first lens.
- **R21**, **R22**: Radii of the second lens.
- **t1**: Thickness of the first lens.
- **t2**: Thickness of the second lens.
- **d**: Spacing between the two lenses.

4. **Update Your Lens Design**: Input these data into your optical system and re-optimize the system.
Note: There is a sample in Sample folder, please view the .docx file to see the performance improvement of this method compared with add a thin plate.
---

## Notes

- This script requires MATLAB with the Global Optimization Toolbox.
- All dimensions are in millimeters unless specified otherwise.
- Ensure you have a working understanding of CODEV and lens design principles.

---

## Contact

For further information or assistance, please feel free to reach out!
