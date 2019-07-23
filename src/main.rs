//Author: Everett Sullivan.
//Date Created: 6-9-2019
//Purpose: expierment with various noise functions and their uses.
//References: https://flafla2.github.io/2014/08/09/perlinnoise.html
//            https://medium.com/100-days-of-algorithms/day-88-perlin-noise-96d23158a44c
//            https://thebookofshaders.com/11/
//            https://www.scratchapixel.com/lessons/procedural-generation-virtual-worlds/perlin-noise-part-2/perlin-noise
//To Do: Figrue out how to seed a random number generator.
//       find out how to make looping 2D perlin noise look better.
//       There seems to be some aterfacting going on, check that alogirthm is implemented correctly
//       Create a method for 2D perlin noise that repeats in two directions.

mod perlin;
use perlin::PerlinGrid;
use perlin::MyNDVec;

mod diamond_square;
use diamond_square::SquareDiamondGrid;

extern crate image;

extern crate clap;
use clap::{Arg, App, SubCommand};

fn main() {

    let matches = App::new("Noise Maker")
        .version("1.0")
        .author("Everett Sullivan")
        .about("Create Perline Noise")
        .arg(Arg::with_name("dimension")
                                    .help("Determines the dimension of the perlin grid")
                                    .index(1)
                                    .required(true))
        .arg(Arg::with_name("image_x_size")
                                    .help("Sets the x dimension of the png file")
                                    .index(2)
                                    .required(true))
        .arg(Arg::with_name("image_y_size")
                                    .help("Sets the y dimension of the png file")
                                    .index(3)
                                    .required(true))
        .arg(Arg::with_name("grid_size")
                                    .help("Sets the size of the perlin grid")
                                    .index(4)
                                    .required(true))
        .arg(Arg::with_name("output_file")
                                    .help("Sets the name of the output file.")
                                    .index(5)
                                    .required(true))
        .arg(Arg::with_name("octaves")
            .help("The number of octaves to use.")
            .takes_value(true)
            .short("o")
            .long("octaves")
            .requires("persistence"))
        .arg(Arg::with_name("persistence")
            .help(".")
            .takes_value(true)
            .short("p")
            .long("persistence")
            .requires("octaves")) //redundant since .requires goes both ways.
        .arg(Arg::with_name("looping")
            .help("Force noise to match up on edges.")
            .takes_value(false)
            .short("l")
            .long("looping"))
        .get_matches();

    //safe to unwrap since the arugment is required.
    let output_file_name: String = matches.value_of("output_file").unwrap().to_string();
    let dimension = matches.value_of("dimension").unwrap().parse::<u32>().unwrap();
    let image_x_size = matches.value_of("image_x_size").unwrap().parse::<u32>().unwrap();
    let image_y_size = matches.value_of("image_y_size").unwrap().parse::<u32>().unwrap();
    let grid_size = matches.value_of("grid_size").unwrap().parse::<usize>().unwrap();
    let octaves;
    let persistence;
    if matches.is_present("octaves") {
        octaves = matches.value_of("octaves").unwrap().parse::<usize>().unwrap();
        persistence = matches.value_of("persistence").unwrap().parse::<f64>().unwrap();
    }else {
        octaves = 1;
        persistence = 1.0;
    }
    let looping = matches.is_present("looping");

    match (dimension,octaves,looping) {
        (_,0,_) => println!("Error"),
        (1,1,false) => create_perlin_noise_1d(image_x_size, image_y_size, grid_size, output_file_name),
        (2,1,false) => create_perlin_noise_2d(image_x_size, image_y_size, grid_size, grid_size, output_file_name),
        (3,1,false) => create_perlin_noise_3d_slice(image_x_size, image_y_size, grid_size, grid_size, output_file_name),
        (1,_,false) => create_perlin_noise_1d_octaves(image_x_size, image_y_size, octaves, persistence,output_file_name),
        (2,_,false) => create_perlin_noise_2d_octaves(image_x_size, image_y_size, octaves, persistence,output_file_name),
        (3,_,false) => create_perlin_noise_3d_slice_octaves(image_x_size, image_y_size, octaves, persistence,output_file_name),
        (1,_,true) => create_perlin_noise_1d_loop(image_x_size, image_y_size, octaves, persistence,output_file_name),
        (2,_,true) => create_perlin_noise_2d_cylinder(image_x_size, image_y_size, octaves, persistence,output_file_name),
        _ => println!("Error"),
    }

    create_diamond_square_noise_2d("squareDiamond.png".to_string())

    //create_perlin_noise_2d(1000,1000,3,5,"test2D_2.png".to_string());

    //create_perlin_noise_2d_octaves(1000,1000,8,0.5,"test2Doct_2.png".to_string());

    //create_perlin_noise_1d(400, 100, 4, output_file_name);

    //create_perlin_noise_1d_octaves(400,100,6,0.7,"test1Doct_2.png".to_string() );

    //create_perlin_noise_1d(400, 100, 4, "blah1.png".to_string());

    //create_perlin_noise_1d_octaves(400,100,1,1.0,"blah2.png".to_string() );

    //create_perlin_noise_3d_slice(1000, 1000, 3, 4, "test3Dslice_2.png".to_string());

    //create_perlin_noise_3d_slice_octaves(1000, 1000, 4, 0.5, "test3Dsliceoct_2.png".to_string());
}

//create_diamond_square_noise_2d
//Purpose:
//    Creates a png of size imgx by imgy image drawing a curve created by 1D perlin noise.
//Pre-conditions:
//    gridx is at least 2, since we need to be able to sample cells.
fn create_diamond_square_noise_2d(file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(1025, 1025);

    let my_noise = SquareDiamondGrid::init_2d(10,0.5,1.0,0.0,0.3,0.7,1.0);
    let mut max_value = 0.0;
    let mut min_value = 2.0;

    for x in 0..1025 {
        for y in 0..1025 {
            let value = my_noise.get_value(vec![x,y]);
            if value > max_value {
                max_value = value;
            }
            if value < min_value {
                min_value = value;
            }
        }
    }

    for x in 0..1025 {
        for y in 0..1025 {
            let value = my_noise.get_value(vec![x,y]);
            let pixel = imgbuf.get_pixel_mut(x as u32, y as u32);
            let mono = (((value - min_value)/(max_value-min_value))*256.0) as u8;
            *pixel = image::Rgb([mono, mono, mono]);
        }
    }

    imgbuf.save(file_name).unwrap();
}

//create_perlin_noise_1d
//Purpose:
//    Creates a png of size imgx by imgy image drawing a curve created by 1D perlin noise.
//Pre-conditions:
//    gridx is at least 2, since we need to be able to sample cells.
fn create_perlin_noise_1d(imgx: u32, imgy: u32, gridx: usize, file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let my_gradiant_grid = PerlinGrid::init(vec![gridx]);

    let grid_to_img_x_ratio = ((gridx-1) as f64)/(imgx as f64);

    for x in 0..imgx {
        let vec = MyNDVec::OneDimVec((0.5 + x as f64)*grid_to_img_x_ratio);
        let raw = my_gradiant_grid.get_perlin(vec);
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

    imgbuf.save(file_name).unwrap();
}

//create_perlin_noise_1d_octaves
//Purpose:
//    Creates a png of size imgx by imgy image drawing a curve created by 1D perlin noise using octaves.
//Pre-conditions:
//    octaves is at least 1.
fn create_perlin_noise_1d_octaves(imgx: u32, imgy: u32, octaves: usize, persistence: f64, file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    let mut grids: Vec<PerlinGrid> = Vec::new();

    let mut double = 1;
    for _i in 0..octaves {
        grids.push(PerlinGrid::init(vec![double+1]));
        double = double*2;
    }

    for x in 0..imgx {
        let mut total_value = 0.0;
        let mut frequency = 1.0;
        let mut amplitude = 1.0;
        let mut max_value = 0.0; //for scaling back to the interval (-1,1).
        for i in 0..octaves {
            let vec = MyNDVec::OneDimVec((0.5 + x as f64)*(frequency)/(imgx as f64));
            total_value = total_value + grids[i].get_perlin(vec) * amplitude;
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

    imgbuf.save(file_name).unwrap();
}

//create_perlin_noise_1d_loop
//Purpose:
//    Creates a png of size imgx by imgy image drawing a curve created by 1D circular silce of a 2D noise
//Pre-conditions:
//    octaves is at least 1.
fn create_perlin_noise_1d_loop(imgx: u32, imgy: u32, octaves: usize, persistence: f64, file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    let mut grids: Vec<PerlinGrid> = Vec::new();

    let mut double = 1;
    for _i in 0..octaves {
        grids.push(PerlinGrid::init(vec![double+1,double+1]));
        double = double*2;
    }

    for x in 0..imgx {
        let mut total_value = 0.0;
        let mut frequency = 1.0;
        let mut amplitude = 1.0;
        let mut max_value = 0.0; //for scaling back to the interval (-1,1).
        for i in 0..octaves {
            let mult_factor = 2.0*std::f64::consts::PI*(0.5 + x as f64)*(frequency)/(imgx as f64);
            let vec = MyNDVec::TwoDimVec(perlin::My2DVec{x: frequency*(1.0 + mult_factor.cos())/2.0, y: frequency*(1.0 + mult_factor.sin())/2.0});
            total_value = total_value + grids[i].get_perlin(vec) * amplitude;
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

    imgbuf.save(file_name).unwrap();
}

//create_perlin_noise_2d
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a 2D perlin noise function for pixels.
//Pre-conditions:
//    gridx and gridy is at least 2, since we need to be able to sample cells.
fn create_perlin_noise_2d(imgx: u32, imgy: u32, gridx: usize, gridy: usize, file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let my_gradiant_grid = PerlinGrid::init(vec![gridx,gridy]);

    let grid_to_img_x_ratio = ((gridx-1) as f64)/(imgx as f64);
    let grid_to_img_y_ratio = ((gridy-1) as f64)/(imgy as f64);

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let vec = MyNDVec::TwoDimVec(perlin::My2DVec{x: (0.5 + x as f64)*grid_to_img_x_ratio, y: (0.5 + y as f64)*grid_to_img_y_ratio});
        let raw = my_gradiant_grid.get_perlin(vec);
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save(file_name).unwrap();
}

//create_perlin_noise_2d_octaves
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a 2D perlin noise function using octaves for pixels.
//Pre-conditions:
//    octaves is at least 1.
fn create_perlin_noise_2d_octaves(imgx: u32, imgy: u32, octaves: usize, persistence: f64, file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    let mut grids: Vec<PerlinGrid> = Vec::new();

    let mut double = 1;
    for _i in 0..octaves {
        grids.push(PerlinGrid::init(vec![double+1,double+1]));
        double = double*2;
    }

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let mut total_value = 0.0;
        let mut frequency = 1.0;
        let mut amplitude = 1.0;
        let mut max_value = 0.0; //for scaling back to the interval (-1,1).
        for i in 0..octaves {
            let vec = MyNDVec::TwoDimVec(perlin::My2DVec{x: (0.5 + x as f64)*(frequency)/(imgx as f64), y: (0.5 + y as f64)*(frequency)/(imgy as f64)});
            total_value = total_value + grids[i].get_perlin(vec) * amplitude;
            max_value = max_value + amplitude;
            amplitude = amplitude * persistence;
            frequency = frequency*2.0;
        }
        let raw = total_value/max_value; // will be a value between -1 and 1 inclusive
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save(file_name).unwrap();
}

//create_perlin_noise_2d_cylinder
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a cylinderical slice of 3D perlin noise.
//Pre-conditions:
//    octaves is at least 1.
fn create_perlin_noise_2d_cylinder(imgx: u32, imgy: u32, octaves: usize, persistence: f64, file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    let mut grids: Vec<PerlinGrid> = Vec::new();

    let mut double = 1;
    for _i in 0..octaves {
        grids.push(PerlinGrid::init(vec![double+1,double+1,double+1]));
        double = double*2;
    }

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let mut total_value = 0.0;
        let mut frequency = 1.0;
        let mut amplitude = 1.0;
        let mut max_value = 0.0; //for scaling back to the interval (-1,1).
        for i in 0..octaves {
            let mult_factor = 2.0*std::f64::consts::PI*(0.5 + x as f64)*(frequency)/(imgx as f64);
            let vec = MyNDVec::ThreeDimVec(perlin::My3DVec{x: frequency*(1.0 + mult_factor.cos())/2.0, y: frequency*(1.0 + mult_factor.cos())/2.0, z: (0.5 + y as f64)*(frequency)/(imgy as f64)});
            total_value = total_value + grids[i].get_perlin(vec) * amplitude;
            max_value = max_value + amplitude;
            amplitude = amplitude * persistence;
            frequency = frequency*2.0;
        }
        let raw = total_value/max_value; // will be a value between -1 and 1 inclusive
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save(file_name).unwrap();
}

//create_perlin_noise_3d_slice
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a 3D perlin noise function for pixels.
//Pre-conditions:
//    gridx and gridy is at least 2, since we need to be able to sample cells.
fn create_perlin_noise_3d_slice(imgx: u32, imgy: u32, gridx: usize, gridy: usize, file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);

    let my_gradiant_grid = PerlinGrid::init(vec![gridx,gridy,2]);

    let grid_to_img_x_ratio = ((gridx-1) as f64)/(imgx as f64);
    let grid_to_img_y_ratio = ((gridy-1) as f64)/(imgy as f64);

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let vec = MyNDVec::ThreeDimVec(perlin::My3DVec{x: (0.5 + x as f64)*grid_to_img_x_ratio, y: (0.5 + y as f64)*grid_to_img_y_ratio, z: 0.5});
        let raw = my_gradiant_grid.get_perlin(vec);
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save(file_name).unwrap();
}

//create_perlin_noise_3d_slice_octaves
//Purpose:
//    Creates a png of size imgx by imgy image by sampleing a 3D perlin noise function using octaves for pixels.
//Pre-conditions:
//    octaves is at least 1.
fn create_perlin_noise_3d_slice_octaves(imgx: u32, imgy: u32, octaves: usize, persistence: f64, file_name: String){
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(imgx, imgy);
    let mut grids: Vec<PerlinGrid> = Vec::new();

    let mut double = 1;
    for _i in 0..octaves {
        grids.push(PerlinGrid::init(vec![double+1,double+1,2]));
        double = double*2;
    }

    // Iterate over the coordinates and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let mut total_value = 0.0;
        let mut frequency = 1.0;
        let mut amplitude = 1.0;
        let mut max_value = 0.0; //for scaling back to the interval (-1,1).
        for i in 0..octaves {
            let vec = MyNDVec::ThreeDimVec(perlin::My3DVec{x: (0.5 + x as f64)*(frequency)/(imgx as f64), y: (0.5 + y as f64)*(frequency)/(imgy as f64), z: 0.5});
            total_value = total_value + grids[i].get_perlin(vec) * amplitude;
            max_value = max_value + amplitude;
            amplitude = amplitude * persistence;
            frequency = frequency*2.0;
        }
        let raw = total_value/max_value;
        let mono = (((raw + 1.0)/2.0)*256.0) as u8;
        *pixel = image::Rgb([mono, mono, mono]);
    }

    imgbuf.save(file_name).unwrap();
}
