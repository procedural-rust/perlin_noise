//Author: Everett Sullivan.
//Date Created: 7-11-2019
//Purpose: Allows creatation and smampling of the square diamond algorithm
//ref http://www.lighthouse3d.com/opengl/terrain/index.php?mpd2
// http://stevelosh.com/blog/2016/06/diamond-square/#edge-cases
//To Do: implement Square-Diamond for 3D

use rand::Rng;
use std::cmp;
use std::f64;

#[derive(Debug, Clone)]
pub struct SquareDiamondGrid1D {
    iterations: usize,
    displacement: f64,
    roughness: f64,
    gradiant_grid: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct SquareDiamondGrid2D {
    iterations: usize,
    displacement: f64,
    roughness: f64,
    gradiant_grid: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct SquareDiamondGrid3D {
    iterations: usize,
    displacement: f64,
    roughness: f64,
    gradiant_grid: Vec<Vec<Vec<f64>>>,
}

#[derive(Debug)]
pub enum SquareDiamondGrid {
    OneDimGrid(SquareDiamondGrid1D),
    TwoDimGrid(SquareDiamondGrid2D),
    ThreeDimGrid(SquareDiamondGrid3D),
}

impl SquareDiamondGrid {

    //init_1d
    //Purpose:
    //    Creates a grid of noise by using the diamond square method (which doesn't have any diamond steps in 1D).
    //    iterations determines the size of the grid, the size will be 2^iterations+1.
    //    displacement will determine how much noise is introduced at each step, modified by 2^-roughness each iteration.
    //    roughness will determine how much variation ther will be in the entire grid, typically 1.0, higher values will be more noisy
    //Pre-conditions:
    //    None
    pub fn init_1d(
        iterations: usize,
        displacement: f64,
        roughness: f64,
        a: f64,
        b: f64,
    ) -> SquareDiamondGrid {
        let size = 2_usize.pow(iterations as u32)+1;
        let mut current_displacement = displacement;
        let mut gradiant_grid = vec![0.0; size];
        gradiant_grid[0] = a;
        gradiant_grid[size-1] = b;
        for i in 0..iterations {
            //square step
            let pow_i = 2_usize.pow(i as u32);
            for x in 0..pow_i {
                //itertate over all bottom left coners of the squares.
                let step = 2_usize.pow((iterations-i) as u32); // distance between squares.
                let my_x = x*step;
                gradiant_grid[my_x + step/2] = (gradiant_grid[my_x]
                                             + gradiant_grid[my_x + step])/2.0
                                             + rand::thread_rng().gen_range(-current_displacement, current_displacement);
            }
            current_displacement *= (2.0_f64).powf(-roughness);

        }
        return SquareDiamondGrid::OneDimGrid(SquareDiamondGrid1D {iterations, displacement, roughness, gradiant_grid });
    }

    //init_2d
    //Purpose:
    //    Creates a grid of noise by using the diamond square method.
    //    iterations determines the size of the grid, the size will be 2^iterations+1.
    //    displacement will determine how much noise is introduced at each step, modified by 2^-roughness each iteration.
    //    roughness will determine how much variation ther will be in the entire grid, typically 1.0, higher values will be more noisy
    //Pre-conditions:
    //    None
    pub fn init_2d(
        iterations: usize,
        displacement: f64,
        roughness: f64,
        aa: f64,
        ab: f64,
        ba: f64,
        bb: f64,
    ) -> SquareDiamondGrid {
        let size = 2_usize.pow(iterations as u32)+1;
        let mut current_displacement = displacement;
        let mut gradiant_grid = vec![vec![0.0; size]; size];
        gradiant_grid[0][0] = aa;
        gradiant_grid[0][size-1] = ab;
        gradiant_grid[size-1][0] = ba;
        gradiant_grid[size-1][size-1] = bb;
        for i in 0..iterations {
            //square step
            let pow_i = 2_usize.pow(i as u32);
            for x in 0..pow_i {
                for y in 0..pow_i {
                    //itertate over all bottom left coners of the squares.
                    let step = 2_usize.pow((iterations-i) as u32); // distance between squares.
                    let my_x = x*step;
                    let my_y = y*step;
                    gradiant_grid[my_x + step/2][my_y + step/2] = (gradiant_grid[my_x][my_y]
                                                                + gradiant_grid[my_x + step][my_y]
                                                                + gradiant_grid[my_x][my_y + step]
                                                                + gradiant_grid[my_x + step][my_y + step])/4.0
                                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                }
            }
            //diamond step
            //diamond steps steps on the edge of the grid only have 3 values to average.
            let half_step = 2_usize.pow((iterations-i-1) as u32); //distance from center of diamond to edge.
            for k in 0..i {
                let my_index = (k*2 + 1)*half_step;
                gradiant_grid[0][my_index] = (gradiant_grid[0][my_index - half_step]
                                           + gradiant_grid[0][my_index + half_step]
                                           + gradiant_grid[half_step][my_index])/3.0
                                           + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[size-1][my_index] = (gradiant_grid[size-1][my_index - half_step]
                                                + gradiant_grid[size-1][my_index + half_step]
                                                + gradiant_grid[size-1-half_step][my_index])/3.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[my_index][0] = (gradiant_grid[my_index - half_step][0]
                                           + gradiant_grid[my_index + half_step][0]
                                           + gradiant_grid[my_index][half_step])/3.0
                                           + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[my_index][size-1] = (gradiant_grid[my_index - half_step][size-1]
                                                + gradiant_grid[my_index + half_step][size-1]
                                                + gradiant_grid[my_index][size-1-half_step])/3.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
            }
            //the rest of the diamond steps have four values to average
            for long_side in 0..pow_i {
                for short_side in 0..(pow_i - 1) {
                    let my_x = long_side*2*half_step + half_step;
                    let my_y = short_side*2*half_step + 2*half_step;
                    gradiant_grid[my_x][my_y] = (gradiant_grid[my_x + half_step][my_y]
                                                + gradiant_grid[my_x - half_step][my_y]
                                                + gradiant_grid[my_x][my_y + half_step]
                                                + gradiant_grid[my_x][my_y - half_step])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                    gradiant_grid[my_y][my_x] = (gradiant_grid[my_y + half_step][my_x]
                                                + gradiant_grid[my_y - half_step][my_x]
                                                + gradiant_grid[my_y][my_x + half_step]
                                                + gradiant_grid[my_y][my_x - half_step])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                }
            }
            current_displacement *= (2.0_f64).powf(-roughness);

        }
        return SquareDiamondGrid::TwoDimGrid(SquareDiamondGrid2D {iterations, displacement, roughness, gradiant_grid });
    }

    //init_3d
    //Purpose:
    //    Creates a grid of noise by using the diamond square method.
    //    iterations determines the size of the grid, the size will be 2^iterations+1.
    //    displacement will determine how much noise is introduced at each step, modified by 2^-roughness each iteration.
    //    roughness will determine how much variation ther will be in the entire grid, typically 1.0, higher values will be more noisy
    //Pre-conditions:
    //    None
    pub fn init_3d(
        iterations: usize,
        displacement: f64,
        roughness: f64,
        aaa: f64,
        aab: f64,
        aba: f64,
        abb: f64,
        baa: f64,
        bab: f64,
        bba: f64,
        bbb: f64,
    ) -> SquareDiamondGrid {
        let size = 2_usize.pow(iterations as u32)+1;
        let mut current_displacement = displacement;
        let mut gradiant_grid = vec![vec![vec![0.0; size]; size]; size];
        gradiant_grid[0][0][0] = aaa;
        gradiant_grid[0][0][size-1] = aab;
        gradiant_grid[0][size-1][0] = aba;
        gradiant_grid[0][size-1][size-1] = abb;
        gradiant_grid[size-1][0][0] = baa;
        gradiant_grid[size-1][0][size-1] = bab;
        gradiant_grid[size-1][size-1][0] = bba;
        gradiant_grid[size-1][size-1][size-1] = bbb;
        for i in 0..iterations {
            //square step
            let pow_i = 2_usize.pow(i as u32);
            for x in 0..pow_i {
                for y in 0..pow_i {
                    for z in 0..pow_i {
                        //itertate over all bottom left coners of the squares.
                        let step = 2_usize.pow((iterations-i) as u32); // distance between squares.
                        let my_x = x*step;
                        let my_y = y*step;
                        let my_z = z*step;
                        gradiant_grid[my_x + step/2][my_y + step/2][my_z + step/2] = (gradiant_grid[my_x][my_y][my_z]
                                                                                   + gradiant_grid[my_x + step][my_y][my_z]
                                                                                   + gradiant_grid[my_x][my_y + step][my_z]
                                                                                   + gradiant_grid[my_x + step][my_y + step][my_z]
                                                                                   + gradiant_grid[my_x][my_y][my_z + step]
                                                                                   + gradiant_grid[my_x + step][my_y][my_z + step]
                                                                                   + gradiant_grid[my_x][my_y + step][my_z + step]
                                                                                   + gradiant_grid[my_x + step][my_y + step][my_z + step])/8.0
                                                                                   + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                    }
                }
            }
            //diamond step
            //diamond steps steps on the edge of the grid only have 4 values to average.
            let half_step = 2_usize.pow((iterations-i-1) as u32); //distance from center of diamond to edge.
            for k in 0..i {
                let my_index = (k*2 + 1)*half_step;
                gradiant_grid[0][0][my_index] = (gradiant_grid[0][0][my_index - half_step]
                                           + gradiant_grid[0][0][my_index + half_step]
                                           + gradiant_grid[0][half_step][my_index]
                                           + gradiant_grid[half_step][0][my_index])/4.0
                                           + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[0][size-1][my_index] = (gradiant_grid[0][size-1][my_index - half_step]
                                                + gradiant_grid[0][size-1][my_index + half_step]
                                                + gradiant_grid[0][size-1-half_step][my_index]
                                                + gradiant_grid[half_step][size-1][my_index])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[size-1][0][my_index] = (gradiant_grid[size-1][0][my_index - half_step]
                                           + gradiant_grid[size-1][0][my_index + half_step]
                                           + gradiant_grid[size-1][half_step][my_index]
                                           + gradiant_grid[size-1-half_step][0][my_index])/4.0
                                           + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[size-1][size-1][my_index] = (gradiant_grid[size-1][size-1][my_index - half_step]
                                                + gradiant_grid[size-1][size-1][my_index + half_step]
                                                + gradiant_grid[size-1][size-1-half_step][my_index]
                                                + gradiant_grid[size-1-half_step][size-1][my_index])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[0][my_index][0] = (gradiant_grid[0][my_index - half_step][0]
                                           + gradiant_grid[0][my_index + half_step][0]
                                           + gradiant_grid[half_step][my_index][0]
                                           + gradiant_grid[0][my_index][half_step])/4.0
                                           + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[0][my_index][size-1] = (gradiant_grid[0][my_index - half_step][size-1]
                                                + gradiant_grid[0][my_index + half_step][size-1]
                                                + gradiant_grid[half_step][my_index][size-1]
                                                + gradiant_grid[0][my_index][size-1-half_step])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[size-1][my_index][0] = (gradiant_grid[size-1][my_index - half_step][0]
                                           + gradiant_grid[size-1][my_index + half_step][0]
                                           + gradiant_grid[size-1-half_step][my_index][0]
                                           + gradiant_grid[size-1][my_index][half_step])/4.0
                                           + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[size-1][my_index][size-1] = (gradiant_grid[size-1][my_index - half_step][size-1]
                                                + gradiant_grid[size-1][my_index + half_step][size-1]
                                                + gradiant_grid[size-1-half_step][my_index][size-1]
                                                + gradiant_grid[size-1][my_index][size-1-half_step])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[my_index][0][0] = (gradiant_grid[my_index - half_step][0][0]
                                           + gradiant_grid[my_index + half_step][0][0]
                                           + gradiant_grid[my_index][0][half_step]
                                           + gradiant_grid[my_index][half_step][0])/4.0
                                           + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[my_index][size-1][0] = (gradiant_grid[my_index - half_step][size-1][0]
                                                + gradiant_grid[my_index + half_step][size-1][0]
                                                + gradiant_grid[my_index][size-1][half_step]
                                                + gradiant_grid[my_index][size-1-half_step][0])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[my_index][0][size-1] = (gradiant_grid[my_index - half_step][0][size-1]
                                           + gradiant_grid[my_index + half_step][0][size-1]
                                           + gradiant_grid[my_index][0][size-1-half_step]
                                           + gradiant_grid[my_index][half_step][size-1])/4.0
                                           + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                gradiant_grid[my_index][size-1][size-1] = (gradiant_grid[my_index - half_step][size-1][size-1]
                                                + gradiant_grid[my_index + half_step][size-1][size-1]
                                                + gradiant_grid[my_index][size-1][size-1-half_step]
                                                + gradiant_grid[my_index][size-1-half_step][size-1])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
            }
            ////////////
            //Not yet done!!!!
            ////////////
            //diamond steps steps on the sides of the grid only have 5 values to average.
            //the rest of the diamond steps have four values to average
            for long_side in 0..pow_i {
                for short_side in 0..(pow_i - 1) {
                    let my_x = long_side*2*half_step + half_step;
                    let my_y = short_side*2*half_step + 2*half_step;
                    gradiant_grid[my_x][my_y][0] = (gradiant_grid[my_x + half_step][my_y][0]
                                                + gradiant_grid[my_x - half_step][my_y][0]
                                                + gradiant_grid[my_x][my_y + half_step][0]
                                                + gradiant_grid[my_x][my_y - half_step][0])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                    gradiant_grid[my_y][my_x][0] = (gradiant_grid[my_y + half_step][my_x][0]
                                                + gradiant_grid[my_y - half_step][my_x][0]
                                                + gradiant_grid[my_y][my_x + half_step][0]
                                                + gradiant_grid[my_y][my_x - half_step][0])/4.0
                                                + rand::thread_rng().gen_range(-current_displacement, current_displacement);
                }
            }
            current_displacement *= (2.0_f64).powf(-roughness);

        }
        return SquareDiamondGrid::ThreeDimGrid(SquareDiamondGrid3D {iterations, displacement, roughness, gradiant_grid });
    }

    pub fn get_value(&self, index: Vec<usize>) -> f64{
        match (self,index.len()) {
            (_,0) => 0.0,
            (SquareDiamondGrid::OneDimGrid(gradiant_grid_1d),_) => gradiant_grid_1d.gradiant_grid[index[0]],
            (SquareDiamondGrid::TwoDimGrid(gradiant_grid_2d),1) => gradiant_grid_2d.gradiant_grid[index[0]][0],
            (SquareDiamondGrid::TwoDimGrid(gradiant_grid_2d),_) => gradiant_grid_2d.gradiant_grid[index[0]][index[1]],
            (SquareDiamondGrid::ThreeDimGrid(gradiant_grid_3d),1) => gradiant_grid_3d.gradiant_grid[index[0]][0][0],
            (SquareDiamondGrid::ThreeDimGrid(gradiant_grid_3d),2) => gradiant_grid_3d.gradiant_grid[index[0]][index[1]][0],
            (SquareDiamondGrid::ThreeDimGrid(gradiant_grid_3d),_) => gradiant_grid_3d.gradiant_grid[index[0]][index[1]][index[2]],
            (_,_) => 0.0,
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diamond_square_test() {

        for i in 0..10 {
            let noise = SquareDiamondGrid1D::init_2d(i,0.00000000000001,0.0,1.0,1.0);
            let size = 2_usize.pow(i as u32)+1;
            for x in 0..size {
                (1.0 - noise.get_value(vec![x])).abs() < 0.01;
            }
        }

        for i in 0..10 {
            let noise = SquareDiamondGrid2D::init_2d(i,0.00000000000001,0.0,1.0,1.0,1.0,1.0);
            let size = 2_usize.pow(i as u32)+1;
            for x in 0..size {
                for y in 0..size {
                    (1.0 - noise.get_value(vec![x,y])).abs() < 0.01;
                }
            }
        }
    }

}