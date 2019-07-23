//Author: Everett Sullivan.
//Date Created: 7-9-2019
//Purpose: Allows creatation and smampling of perlin noise algorithms.
//To Do: Implment Perline 4D sampling
//       Better Error checking.

use rand::Rng;
use std::cmp;
use std::f64;

#[derive(Debug, Copy, Clone,PartialEq,PartialOrd)]
pub struct My2DVec {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Copy, Clone,PartialEq,PartialOrd)]
pub struct My3DVec {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone)]
pub struct GradiantGrid1D {
    width: usize,
    gradiant_grid: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct GradiantGrid2D {
    width: usize,
    height: usize,
    gradiant_grid: Vec<Vec<My2DVec>>,
}

#[derive(Debug, Clone)]
pub struct GradiantGrid3D {
    width: usize,
    height: usize,
    depth: usize,
    gradiant_grid: Vec<Vec<Vec<My3DVec>>>,
}

#[derive(Debug)]
pub enum PerlinGrid {
    OneDimGrid(GradiantGrid1D),
    TwoDimGrid(GradiantGrid2D),
    ThreeDimGrid(GradiantGrid3D),
}

#[derive(Debug,Copy,Clone)]
pub enum MyNDVec  {
    OneDimVec(f64),
    TwoDimVec(My2DVec),
    ThreeDimVec(My3DVec),
}

impl PerlinGrid {

    //init
    //Purpose:
    //    Creates a Perline grid of the given dimension.
    //    The length of the vector determines the dimension of the grid.
    //    If the length is zero, a default grid of dimension one with size 2 is created.
    //    If the length is four or more the extra entries are ingored.
    //    If any of the entries in the vector is less than 2, they are replaced with 2.
    //Pre-conditions:
    //    the length of the vector is between 1 and 3.
    //    Each entry in the vector must be at least two since otherwise there will be no cells for sampling.
    pub fn init(
        grid_dimensions: Vec<usize>,
    ) -> PerlinGrid {
        match grid_dimensions.len() {
            0 => PerlinGrid::OneDimGrid(create_one_dim_grid(2,random_1d_gradiant_vector_1)),
            1 => PerlinGrid::OneDimGrid(create_one_dim_grid(cmp::max(grid_dimensions[0], 2),random_1d_gradiant_vector_1)),
            2 => PerlinGrid::TwoDimGrid(create_two_dim_grid(cmp::max(grid_dimensions[0], 2),cmp::max(grid_dimensions[1], 2),random_2d_gradiant_vector_3)),
            _ => PerlinGrid::ThreeDimGrid(create_three_dim_grid(cmp::max(grid_dimensions[0], 2),cmp::max(grid_dimensions[1], 2),cmp::max(grid_dimensions[2], 2),random_3d_gradiant_vector_1)),
        }
    }

    //get_perlin
    //Purpose:
    //    Returns the noise value at vector of floats using Perlin noise.
    //    If the dimension of the vector and Perlin grid are not the same, the method returns zero.
    //Pre-conditions:
    //    The dimension of the vector and Perlin grid must be the same.
    pub fn get_perlin(&self, index: MyNDVec) -> f64{
        match (self,index) {
            (PerlinGrid::OneDimGrid(gradiant_grid_1d),MyNDVec::OneDimVec(my_1d_vec)) => get_one_dim_perlin(&gradiant_grid_1d, my_1d_vec),
            (PerlinGrid::TwoDimGrid(gradiant_grid_2d),MyNDVec::TwoDimVec(my_2d_vec)) => get_two_dim_perlin(&gradiant_grid_2d, my_2d_vec.x, my_2d_vec.y),
            (PerlinGrid::ThreeDimGrid(gradiant_grid_3d),MyNDVec::ThreeDimVec(my_3d_vec)) => get_three_dim_perlin(&gradiant_grid_3d, my_3d_vec.x, my_3d_vec.y, my_3d_vec.z),
            (_,_) => 0.0,
        }
    }

}

////////////////////
//Grid creation
////////////////////

//create_one_dim_grid
//Purpose:
//    Returns a 1D array of 1D vectors (f64 values) for Perlin noise.
//    The vectors are randomly chosen using random_gradiant_function.
//Pre-conditions:
//    The variable width is at least 2.
//Notes:
//    The width must be at least two since otherwise there will be no cells from which to sample.
fn create_one_dim_grid(
    width: usize,
    random_gradiant_function: fn() -> f64
) -> GradiantGrid1D {
    let mut gradiant_grid = vec![0.0; width];
    for i in 0..width {
        gradiant_grid[i] = random_gradiant_function();
    }
    return GradiantGrid1D{ width, gradiant_grid};
}

//create_two_dim_grid
//Purpose:
//    Returns a 2D array of 2D vectors for Perlin noise.
//    The vectors are randomly chosen using random_gradiant_function.
//Pre-conditions:
//    The variables width and height are at least 2.
//Notes:
//    The width and height must be at least two since otherwise there will be no cells from which to sample.
fn create_two_dim_grid(
    width: usize,
    height: usize,
    random_gradiant_function: fn() -> My2DVec
) -> GradiantGrid2D {
    let mut gradiant_grid = vec![vec![My2DVec{ x: 0.0, y: 0.0,}; height]; width];
    for row in 0..width {
        for col in 0..height {
            gradiant_grid[row][col] = random_gradiant_function();
        }
    }
    return GradiantGrid2D{ width, height, gradiant_grid};
}

//create_three_dim_grid
//Purpose:
//    Returns all a 3D array of 3D vectors for Perlin noise.
//    The vectors are randomly chosen using random_gradiant_function.
//Pre-conditions:
//    The variables width, height, and depth are at least 2.
//Notes:
//    The width, height, and depth must be at least two since otherwise there will be no cells from which to sample.
fn create_three_dim_grid(
    height: usize,
    width: usize,
    depth: usize,
    random_gradiant_function: fn() -> My3DVec
) -> GradiantGrid3D {
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

////////////////////
//Perlin sampling
////////////////////

//get_one_dim_perlin
//Purpose:
//    Returns the noise value at x using 1D advance perlin noise algorithm.
//Pre-conditions:
//    None
fn get_one_dim_perlin(grid: &GradiantGrid1D, x: f64) -> f64{
    let my_x = x % ((grid.width-1) as f64);// if the point x is outside the gradiant grid we use modulo to fit it inside the grid.
    let cell_x = my_x.floor() as usize; // get integer c where c <= my_x < c + 1
    let mut px = my_x - my_x.floor();
    let a = grid.gradiant_grid[cell_x];
    let b = grid.gradiant_grid[cell_x + 1];
    let a_dot = px*a;
    let b_dot = (1.0-px)*b;
    px = perlin_fade(px);
    let value = linear_interpolation(a_dot,b_dot,px);
    return value;
}

//get_two_dim_perlin
//Purpose:
//    Returns the noise value at (x,y) using the advance perlin noise algorithm.
//Pre-conditions:
//    None
fn get_two_dim_perlin(grid: &GradiantGrid2D, x: f64, y: f64) -> f64{
    let my_x = x % ((grid.width-1) as f64);
    let my_y = y % ((grid.height-1) as f64); // if the point (x,y) is outside the gradiant grid we use modulo to fit it inside the grid.
    let cell_x = my_x.floor() as usize;
    let cell_y = my_y.floor() as usize; // get (c_x,c_y) where c_x <= x < c_x + 1 and c_y <= y < c_y.
    let mut px = my_x - my_x.floor();
    let mut py = my_y - my_y.floor();
    // get gradiants at the four corners of the cell
    let aa = grid.gradiant_grid[cell_x][cell_y];
    let ab = grid.gradiant_grid[cell_x][cell_y + 1];
    let ba = grid.gradiant_grid[cell_x + 1][cell_y];
    let bb = grid.gradiant_grid[cell_x + 1][cell_y + 1];
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

//get_three_dim_perlin
//Purpose:
//    Returns the noise value at (x,y,z) using the advance perlin noise algorithm.
//Pre-conditions:
//    None
fn get_three_dim_perlin(grid: &GradiantGrid3D, x: f64, y: f64, z: f64) -> f64{
    let my_x = x % ((grid.height-1) as f64);
    let my_y = y % ((grid.width-1) as f64);
    let my_z = z % ((grid.depth-1) as f64); // if the point (x,y,z) is outside the gradiant grid we use modulo to fit it inside the grid.
    let cell_x = my_x.floor() as usize;
    let cell_y = my_y.floor() as usize;
    let cell_z = my_z.floor() as usize;
    let mut px = my_x - my_x.floor();
    let mut py = my_y - my_y.floor();
    let mut pz = my_z - my_z.floor();
    // get gradiants at the eight corners of the cell
    let aaa = grid.gradiant_grid[cell_x    ][cell_y    ][cell_z    ];
    let aab = grid.gradiant_grid[cell_x    ][cell_y    ][cell_z + 1];
    let aba = grid.gradiant_grid[cell_x    ][cell_y + 1][cell_z    ];
    let baa = grid.gradiant_grid[cell_x + 1][cell_y    ][cell_z    ];
    let bba = grid.gradiant_grid[cell_x + 1][cell_y + 1][cell_z    ];
    let bab = grid.gradiant_grid[cell_x + 1][cell_y    ][cell_z + 1];
    let abb = grid.gradiant_grid[cell_x    ][cell_y + 1][cell_z + 1];
    let bbb = grid.gradiant_grid[cell_x + 1][cell_y + 1][cell_z + 1];
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

////////////////////
//Random vectoring sampling
////////////////////

//random_1d_gradiant_vector_1
//Purpose:
//    Returns a 1D vector uniformly at random from the interval (-2,2).
//Pre-conditions:
//    None
fn random_1d_gradiant_vector_1() -> f64{
    rand::thread_rng().gen_range(-2.0, 2.0)
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

////////////////////
//Interpolation Functions
////////////////////

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