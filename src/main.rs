//Author: Everett Sullivan.
//Date Created: 6-9-2019
//Purpose: expierment with various noise functions and their uses.
//References: https://flafla2.github.io/2014/08/09/perlinnoise.html
//            https://medium.com/100-days-of-algorithms/day-88-perlin-noise-96d23158a44c
//            https://thebookofshaders.com/11/
//            https://www.scratchapixel.com/lessons/procedural-generation-virtual-worlds/perlin-noise-part-2/perlin-noise
//To Do: Figrue out how to seed a random number generator.
extern crate image;

use rand::Rng;
use std::f64;

#[derive(Debug, Clone)]
struct GradiantGrid1D {
    width: usize,
    gradiant_grid: Vec<f64>,
}

impl GradiantGrid1D{

    //init
    //Purpose:
    //    Returns all a 1D array of 1D vectors (f64 values) for Perlin noise.
    //    The vectors are randomly chosen using random_gradiant_function.
    //Pre-conditions:
    //    The variable width is at least 2.
    //Notes:
    //    The width must be at least two since otherwise there will be no cells from which to sample.
    pub fn init(
        width: usize,
        random_gradiant_function: fn() -> f64
    ) -> GradiantGrid1D {
        if width < 2 {
            panic!("Attempted creation of a grid with no cells!");
        }
        let mut gradiant_grid = vec![0.0; width];
        for i in 0..width {
            gradiant_grid[i] = random_gradiant_function();
        }
        return GradiantGrid1D{ width, gradiant_grid};
    }

    //get_gradiant
    //Purpose:
    //    Returns the gradiant vector at gradiant_grid[x].
    //Pre-conditions:
    //    The variable x are within bounds.
    pub fn get_gradiant(&self, i: usize) -> f64{
		match self.gradiant_grid.get(i) {
			Some(entry) => {return *entry},
			None => {panic!("Attempted access of a grid out of bounds!");},
		}
	}

    //get_perlin
    //Purpose:
    //    Returns the noise value at x using 1D advance perlin noise algorithm.
    //Pre-conditions:
    //    None
    fn get_perlin(&self, x: f64) -> f64{
        let my_x = x % ((self.width-1) as f64);// if the point (x,y) is outside the gradiant grid we use modulo to fit it inside the grid.
        let cell_x = my_x.floor() as usize; // get integer c where c <= my_x < c + 1
        let mut px = my_x - my_x.floor();
        let a = self.get_gradiant(cell_x);
        let b = self.get_gradiant(cell_x + 1);
        let a_dot = px*a;
        let b_dot = (1.0-px)*b;
        px = perlin_fade(px);
        let value = linear_interpolation(a_dot,b_dot,px);
        return value;
    }

}

//random_1d_gradiant_vector_1
//Purpose:
//    Returns a 1D vector uniformly at random from the interval (-2,2).
//Pre-conditions:
//    None
fn random_1d_gradiant_vector_1() -> f64{
    rand::thread_rng().gen_range(-2.0, 2.0)
}

#[derive(Debug, Copy, Clone,PartialEq,PartialOrd)]
struct My2DVec {
    x: f64,
    y: f64,
}

#[derive(Debug, Clone)]
struct GradiantGrid2D {
    width: usize,
    height: usize,
    gradiant_grid: Vec<Vec<My2DVec>>,
}

impl GradiantGrid2D{

    //init
    //Purpose:
    //    Returns all a 2D array of 2D vectors for Perlin noise.
    //    The vectors are randomly chosen using random_gradiant_function.
    //Pre-conditions:
    //    The variables width and height are at least 2.
    //Notes:
    //    The width and height must be at least two since otherwise there will be no cells from which to sample.
    pub fn init(
        width: usize,
        height: usize,
        random_gradiant_function: fn() -> My2DVec
    ) -> GradiantGrid2D {
        if (width < 2) || (height < 2) {
            panic!("Attempted creation of a grid with no height or no width!");
        }
        let mut gradiant_grid = vec![vec![My2DVec{ x: 0.0, y: 0.0,}; width]; height];
        for row in 0..height {
            for col in 0..width {
                gradiant_grid[row][col] = random_gradiant_function();
            }
        }
        return GradiantGrid2D{ width, height, gradiant_grid};
    }

    //get_gradiant
    //Purpose:
    //    Returns the gradiant vector at gradiant_grid[x][y].
    //Pre-conditions:
    //    The variables x and y are within bounds.
    pub fn get_gradiant(&self, x: usize, y:usize) -> My2DVec{
	    match self.gradiant_grid.get(y) {
		    Some(row) => {
			    match row.get(x) {
				    Some(column) => {return *column},
					None => {panic!("Attempted access of a grid out of bounds!");},
				}
			}
			None => {panic!("Attempted access of a grid out of bounds!");},
		}
	}

    //get_perlin
    //Purpose:
    //    Returns the noise value at (x,y) using the advance perlin noise algorithm.
    //Pre-conditions:
    //    None
    fn get_perlin(&self, x: f64, y: f64) -> f64{
        let my_x = x % ((self.width-1) as f64);
        let my_y = y % ((self.height-1) as f64); // if the point (x,y) is outside the gradiant grid we use modulo to fit it inside the grid.
        let cell_x = my_x.floor() as usize;
        let cell_y = my_y.floor() as usize; // get (c_x,c_y) where c_x <= x < c_x + 1 and c_y <= y < c_y.
        let mut px = my_x - my_x.floor();
        let mut py = my_y - my_y.floor();
        // get gradiants at the four corners of the cell
        let aa = self.get_gradiant(cell_x,cell_y);
        let ab = self.get_gradiant(cell_x,cell_y + 1);
        let ba = self.get_gradiant(cell_x + 1,cell_y);
        let bb = self.get_gradiant(cell_x +1 ,cell_y + 1);
        // compute dot products of each gradiant vector with the vector from it's corner to the point.
        let aa_dot = px*aa.x + py*aa.y;
        let ab_dot = px*ab.x + (1.0-py)*ab.y;
        let ba_dot = (1.0-px)*ba.x + py*ba.y;
        let bb_dot = (1.0-px)*bb.x + (1.0-py)*bb.y;
        px = perlin_fade(px);
        py = perlin_fade(py);
        let x1 = linear_interpolation(aa_dot,ba_dot,px);
        let x2 = linear_interpolation(ab_dot,bb_dot,px);
        let value = linear_interpolation(x1,x2,py);
        return value;
    }

}

//random_2d_gradiant_vector_1
//Purpose:
//    Returns a 2D vector uniformly at random from (1,1),(1,-1),(-1,1), and (-1,-1).
//Pre-conditions:
//    None
fn random_2d_gradiant_vector_1() -> My2DVec{
    let vector_choice = rand::thread_rng().gen_range(0, 5);
    match vector_choice {  
        0 => My2DVec{ x:  1.0, y:  1.0,},
        1 => My2DVec{ x:  1.0, y: -1.0,},
        2 => My2DVec{ x: -1.0, y:  1.0,},
        _ => My2DVec{ x: -1.0, y: -1.0,},
    }
}

//random_2d_gradiant_vector_2
//Purpose:
//    Returns a 2D vector uniformly at random from (1,0),(-1,-0),(0,1), and (0,-1).
//Pre-conditions:
//    None
fn random_2d_gradiant_vector_2() -> My2DVec{
    let vector_choice = rand::thread_rng().gen_range(0, 5);
    match vector_choice {  
        0 => My2DVec{ x: 1.0,  y:  0.0,},
        1 => My2DVec{ x: -1.0, y: 0.0,},
        2 => My2DVec{ x: 0.0,  y:  1.0,},
        _ => My2DVec{ x: 0.0,  y: -1.0,},
    }
}

//random_2d_gradiant_vector_3
//Purpose:
//    Returns a 2D vector uniformly at random from the unit circle.
//Pre-conditions:
//    None
fn random_2d_gradiant_vector_3() -> My2DVec{
    let vector_choice = rand::thread_rng().gen_range(0.0, 2.0*f64::consts::PI);
    My2DVec{ x: vector_choice.cos(),  y: vector_choice.sin(),}
}

#[derive(Debug, Copy, Clone,PartialEq,PartialOrd)]
struct My3DVec {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Debug, Clone)]
struct GradiantGrid3D {
    width: usize,
    height: usize,
    depth: usize,
    gradiant_grid: Vec<Vec<Vec<My3DVec>>>,
}

impl GradiantGrid3D{

    //init
    //Purpose:
    //    Returns all a 3D array of 3D vectors for Perlin noise.
    //    The vectors are randomly chosen using random_gradiant_function.
    //Pre-conditions:
    //    The variables width, height, and depth are at least 2.
    //Notes:
    //    The width, height, and depth must be at least two since otherwise there will be no cells from which to sample.
    pub fn init(
        height: usize,
        width: usize,
        depth: usize,
        random_gradiant_function: fn() -> My3DVec
    ) -> GradiantGrid3D {
        if (width < 2) || (height < 2) || (depth < 2) {
            panic!("Attempted creation of a grid with no height or no width!");
        }
        let mut gradiant_grid = vec![vec![vec![My3DVec{ x: 0.0, y: 0.0, z: 0.0}; depth]; width]; height];
        for row in 0..height {
            for col in 0..width {
                for dep in 0..depth {
                    gradiant_grid[row][col][dep] = random_gradiant_function();
                }
            }
        }
        return GradiantGrid3D{ width, height, depth, gradiant_grid};
    }

    //get_gradiant
    //Purpose:
    //    Returns the gradiant vector at gradiant_grid[x][y][z].
    //Pre-conditions:
    //    The variables x,y and z are within bounds.
    pub fn get_gradiant(&self, x: usize, y:usize, z:usize) -> My3DVec{
 	    match self.gradiant_grid.get(x) {
		    Some(row) => {
			    match row.get(y) {
				    Some(column) => {
                        match column.get(z) {
                            Some(entry) => {return *entry},
                            None => {panic!("Attempted access of a grid out of bounds!");},
                        }
                    }
					None => {panic!("Attempted access of a grid out of bounds!");},
				}
			}
			None => {panic!("Attempted access of a grid out of bounds!");},
		}
	}

    //get_perlin
    //Purpose:
    //    Returns the noise value at (x,y,z) using the advance perlin noise algorithm.
    //Pre-conditions:
    //    None
    fn get_perlin(&self, x: f64, y: f64, z: f64) -> f64{
        let my_x = x % ((self.height-1) as f64);
        let my_y = y % ((self.width-1) as f64);
        let my_z = z % ((self.depth-1) as f64); // if the point (x,y,z) is outside the gradiant grid we use modulo to fit it inside the grid.
        let cell_x = my_x.floor() as usize;
        let cell_y = my_y.floor() as usize;
        let cell_z = my_z.floor() as usize;
        let mut px = my_x - my_x.floor();
        let mut py = my_y - my_y.floor();
        let mut pz = my_z - my_z.floor();
        // get gradiants at the eight corners of the cell
        let aaa = self.get_gradiant(cell_x    ,cell_y    ,cell_z);
        let aab = self.get_gradiant(cell_x    ,cell_y    ,cell_z + 1);
        let aba = self.get_gradiant(cell_x    ,cell_y + 1,cell_z);
        let baa = self.get_gradiant(cell_x + 1,cell_y    ,cell_z);
        let bba = self.get_gradiant(cell_x + 1,cell_y + 1,cell_z);
        let bab = self.get_gradiant(cell_x + 1,cell_y    ,cell_z + 1);
        let abb = self.get_gradiant(cell_x    ,cell_y + 1,cell_z + 1);
        let bbb = self.get_gradiant(cell_x + 1,cell_y + 1,cell_z + 1);
        // compute dot products of each gradiant vector with the vector from it's corner to the point.
        let aaa_dot = px*aaa.x       + py*aaa.y       + pz*aaa.z;
        let aab_dot = px*aab.x       + py*aab.y       + (1.0-pz)*aab.z;
        let aba_dot = px*aba.x       + (1.0-py)*aba.y + pz*aba.z;
        let baa_dot = (1.0-px)*baa.x + py*baa.y       + pz*baa.z;
        let bba_dot = (1.0-px)*bba.x + (1.0-py)*bba.y + pz*bba.z;
        let bab_dot = (1.0-px)*bab.x + py*bab.y       + (1.0-pz)*bab.z;
        let abb_dot = px*abb.x       + (1.0-py)*abb.y + (1.0-pz)*abb.z;
        let bbb_dot = (1.0-px)*bbb.x + (1.0-py)*bbb.y + (1.0-pz)*bbb.z;
        px = perlin_fade(px);
        py = perlin_fade(py);
        pz = perlin_fade(pz);
        let x1 = linear_interpolation(aaa_dot,baa_dot,px);
        let x2 = linear_interpolation(aba_dot,bba_dot,px);
        let x3 = linear_interpolation(aab_dot,bab_dot,px);
        let x4 = linear_interpolation(abb_dot,bbb_dot,px);
        let y1 = linear_interpolation(x1,x2,py);
        let y2 = linear_interpolation(x3,x4,py);
        let value = linear_interpolation(y1,y2,pz);
        return value;
    }

}

//random_3d_gradiant_vector_1
//Purpose:
//    Returns a 3D vector uniformly at random from (1,1,0),(−1,1,0),(1,−1,0),(−1,−1,0),(1,0,1),(−1,0,1),(1,0,−1),(−1,0,−1),(0,1,1),(0,−1,1),(0,1,−1),(0,−1,−1)
//Pre-conditions:
//    None
fn random_3d_gradiant_vector_1() -> My3DVec{
    let vector_choice = rand::thread_rng().gen_range(0, 12);
    match vector_choice {  
        0 => My3DVec{ x:  1.0, y:  1.0, z:  0.0},
        1 => My3DVec{ x:  -1.0, y:  1.0, z:  0.0},
        2 => My3DVec{ x:  1.0, y:  -1.0, z:  0.0},
        4 => My3DVec{ x:  -1.0, y:  -1.0, z:  0.0},
        5 => My3DVec{ x:  1.0, y:  0.0, z:  1.0},
        6 => My3DVec{ x:  -1.0, y:  0.0, z:  1.0},
        7 => My3DVec{ x:  1.0, y:  0.0, z:  -1.0},
        8 => My3DVec{ x:  -1.0, y:  0.0, z:  -1.0},
        9 => My3DVec{ x:  0.0, y:  1.0, z:  1.0},
        10 => My3DVec{ x:  0.0, y:  -1.0, z:  1.0},
        11 => My3DVec{ x:  0.0, y:  1.0, z:  -1.0},
        _ => My3DVec{ x:  0.0, y:  -1.0, z:  -1.0},
    }
}

//perlin_fade
//Purpose:
//    Implements the fade function for Perlin's noise
//Pre-conditions:
//    The variable t is between 0 and 1 inclusive.
//fn perlin_fade(t: f64) -> f64{
//   t * t * (3.0 - 2.0 * t) // 3t^2 - 2t^2
//}

//perlin_fade
//Purpose:
//    Implements the fade function for Perlin's advance noise
//Pre-conditions:
//    The variable t is between 0 and 1 inclusive.
fn perlin_fade(t: f64) -> f64{
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0) // 6t^5 - 15t^4 + 10t^3
}

//linear_interpolation
//Purpose:
//    Implements linear interpolation between two floats.
//Pre-conditions:
//    The variable t is between 0 and 1 inclusive.
fn linear_interpolation(x1: f64, x2: f64, t: f64) -> f64{
    x1*(1.0-t) + x2*t
}

//create_perlin_noise_3d_slice
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a 3D perlin noise function for pixels.
//Pre-conditions:
//    gridx and gridy is at least 2, since we need to be able to sample cells.
fn create_perlin_noise_3d_slice(imgx: u32, imgy: u32, gridx: usize, gridy: usize,random_gradiant_function: fn() -> My3DVec){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let my_gradiant_grid = GradiantGrid3D::init(gridx,gridy,2,random_gradiant_function);

    let grid_to_img_x_ratio = ((gridx-1) as f64)/(imgx as f64);
    let grid_to_img_y_ratio = ((gridy-1) as f64)/(imgy as f64);

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let raw = my_gradiant_grid.get_perlin((0.5 + x as f64)*grid_to_img_x_ratio, (0.5 + y as f64)*grid_to_img_y_ratio, 0.5);
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save("test3Dslice.png").unwrap();
}

//create_perlin_noise_3d_slice_octaves
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a 3D perlin noise function using octaves for pixels.
//Pre-conditions:
//    octaves is at least 1.
fn create_perlin_noise_3d_slice_octaves(imgx: u32, imgy: u32, octaves: usize, persistence: f64,random_gradiant_function: fn() -> My3DVec){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    let mut grids: Vec<GradiantGrid3D> = Vec::new();

    let mut double = 1;
    for _i in 0..octaves {
        grids.push(GradiantGrid3D::init(double+1,double+1,2,random_gradiant_function));
        double = double*2;
    }

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let mut total_value = 0.0;
        let mut frequency = 1.0;
        let mut amplitude = 1.0;
        let mut max_value = 0.0; //for scaling back to the interval (-1,1).
        for i in 0..octaves {
            total_value = total_value + grids[i].get_perlin((0.5 + x as f64)*(frequency)/(imgx as f64),(0.5 + y as f64)*(frequency)/(imgy as f64), 0.5) * amplitude;
            max_value = max_value + amplitude;
            amplitude = amplitude * persistence;
            frequency = frequency*2.0;
        }
        let raw = total_value/max_value;
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save("test3Dsliceoct.png").unwrap();
}

//create_perlin_noise_2d
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a 2D perlin noise function for pixels.
//Pre-conditions:
//    gridx and gridy is at least 2, since we need to be able to sample cells.
fn create_perlin_noise_2d(imgx: u32, imgy: u32, gridx: usize, gridy: usize, random_gradiant_function: fn() -> My2DVec){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let my_gradiant_grid = GradiantGrid2D::init(gridx,gridy,random_gradiant_function);

    let grid_to_img_x_ratio = ((gridx-1) as f64)/(imgx as f64);
    let grid_to_img_y_ratio = ((gridy-1) as f64)/(imgy as f64);

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let raw = my_gradiant_grid.get_perlin((0.5 + x as f64)*grid_to_img_x_ratio, (0.5 + y as f64)*grid_to_img_y_ratio);
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save("test2D.png").unwrap();
}

//create_perlin_noise_2d_octaves
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a 2D perlin noise function using octaves for pixels.
//Pre-conditions:
//    octaves is at least 1.
fn create_perlin_noise_2d_octaves(imgx: u32, imgy: u32, octaves: usize, persistence: f64, random_gradiant_function: fn() -> My2DVec){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    let mut grids: Vec<GradiantGrid2D> = Vec::new();

    let mut double = 1;
    for _i in 0..octaves {
        grids.push(GradiantGrid2D::init(double+1,double+1,random_gradiant_function));
        double = double*2;
    }

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let mut total_value = 0.0;
        let mut frequency = 1.0;
        let mut amplitude = 1.0;
        let mut max_value = 0.0; //for scaling back to the interval (-1,1).
        for i in 0..octaves {
            total_value = total_value + grids[i].get_perlin((0.5 + x as f64)*(frequency)/(imgx as f64),(0.5 + y as f64)*(frequency)/(imgy as f64)) * amplitude;
            max_value = max_value + amplitude;
            amplitude = amplitude * persistence;
            frequency = frequency*2.0;
        }
        let raw = total_value/max_value;
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save("test2Doct.png").unwrap();
}

//create_perlin_noise_1d
//Purpose:
//    Creates a png of size imgx by imgy image drawing a curve created by 1D perlin noise.
//Pre-conditions:
//    gridx is at least 2, since we need to be able to sample cells.
fn create_perlin_noise_1d(imgx: u32, imgy: u32, gridx: usize,random_gradiant_function: fn() -> f64){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let my_gradiant_grid = GradiantGrid1D::init(gridx,random_gradiant_function);

    let grid_to_img_x_ratio = ((gridx-1) as f64)/(imgx as f64);

    for x in 0..imgx {
        let raw = my_gradiant_grid.get_perlin((0.5 + x as f64)*grid_to_img_x_ratio);
        let height = (raw + 1.0)*0.5*(imgy as f64);
        for y in 0..imgy {
            let pixel = imgbuf.get_pixel_mut(x, y);
            if (height - (y as f64)).abs() < 1.0 {
                *pixel = image::Rgb([255,255,255]);
            }else{
                *pixel = image::Rgb([0,0,0]);
            }
        }
    }

    imgbuf.save("test1D.png").unwrap();
}

//create_perlin_noise_1d_octaves
//Purpose:
//    Creates a png of size imgx by imgy image drawing a curve created by 1D perlin noise using octaves.
//Pre-conditions:
//    octaves is at least 1.
fn create_perlin_noise_1d_octaves(imgx: u32, imgy: u32, octaves: usize, persistence: f64,random_gradiant_function: fn() -> f64){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    let mut grids: Vec<GradiantGrid1D> = Vec::new();

    let mut double = 1;
    for _i in 0..octaves {
        grids.push(GradiantGrid1D::init(double+1,random_gradiant_function));
        double = double*2;
    }

    for x in 0..imgx {
        let mut total_value = 0.0;
        let mut frequency = 1.0;
        let mut amplitude = 1.0;
        let mut max_value = 0.0; //for scaling back to the interval (-1,1).
        for i in 0..octaves {
            total_value = total_value + grids[i].get_perlin((0.5 + x as f64)*(frequency)/(imgx as f64)) * amplitude;
            max_value = max_value + amplitude;
            amplitude = amplitude * persistence;
            frequency = frequency*2.0;
        }
        let height = (total_value/max_value + 1.0)*0.5*(imgy as f64);
        for y in 0..imgy {
            let pixel = imgbuf.get_pixel_mut(x, y);
            if (height - (y as f64)).abs() < 1.0 {
                *pixel = image::Rgb([255,255,255]);
            }else{
                *pixel = image::Rgb([0,0,0]);
            }
        }
    }

    imgbuf.save("test1Doct.png").unwrap();
}

fn main() {
    //create_perlin_noise_2d(1000,1000,3,5,random_2d_gradiant_vector_3);

    //create_perlin_noise_2d_octaves(1000,1000,8,0.5,random_2d_gradiant_vector_3);

    //create_perlin_noise_1d(400, 100, 4,random_1d_gradiant_vector_1);

    //create_perlin_noise_1d_octaves(400,100,6,0.7,random_1d_gradiant_vector_1);

    create_perlin_noise_3d_slice(1000, 1000, 3, 4, random_3d_gradiant_vector_1);

    create_perlin_noise_3d_slice_octaves(1000, 1000, 4, 0.5, random_3d_gradiant_vector_1);
}
