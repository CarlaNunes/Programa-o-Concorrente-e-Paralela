#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

double c_x_min;
double c_x_max;
double c_y_min;
double c_y_max;

double pixel_width;
double pixel_height;

int iteration_max = 200;

int image_size;
// Array para receber buffers de cada thread
unsigned char ***image_collection;

int i_x_max;
int i_y_max;
int image_buffer_size;

int gradient_size = 16;
int colors[17][3] = {
                        {66, 30, 15},
                        {25, 7, 26},
                        {9, 1, 47},
                        {4, 4, 73},
                        {0, 7, 100},
                        {12, 44, 138},
                        {24, 82, 177},
                        {57, 125, 209},
                        {134, 181, 229},
                        {211, 236, 248},
                        {241, 233, 191},
                        {248, 201, 95},
                        {255, 170, 0},
                        {204, 128, 0},
                        {153, 87, 0},
                        {106, 52, 3},
                        {16, 16, 16},
                    };

// Numero de threads como argv[6] para iterarmos via bash
int nthreads = 1;

void allocate_image_collection(){
    image_collection = (unsigned char ***) malloc(sizeof(unsigned char ***) * nthreads);
};

void init(int argc, char *argv[]){
    if(argc < 6){
        printf("usage: ./mandelbrot_omp c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_omp -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:      ./mandelbrot_omp -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:      ./mandelbrot_omp 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: ./mandelbrot_omp -0.188 -0.012 0.554 0.754 11500\n");
        exit(0);
    }
    else{
        sscanf(argv[1], "%lf", &c_x_min);
        sscanf(argv[2], "%lf", &c_x_max);
        sscanf(argv[3], "%lf", &c_y_min);
        sscanf(argv[4], "%lf", &c_y_max);
        sscanf(argv[5], "%d", &image_size);

        i_x_max           = image_size;
        i_y_max           = image_size;
        image_buffer_size = image_size * image_size;

        pixel_width       = (c_x_max - c_x_min) / i_x_max;
        pixel_height      = (c_y_max - c_y_min) / i_y_max;
    };
};

void write_to_file(){
    FILE * file;
    char * filename               = "output.ppm";
    char * comment                = "# ";

    int max_color_component_value = 255;

    file = fopen(filename,"wb");

    fprintf(file, "P6\n %s\n %d\n %d\n %d\n", comment,
            i_x_max, i_y_max, max_color_component_value);

    // Variavel que coordena a divisão do tamanho dos buffers e divisao de trabalho.
    int buffer_division;
    if(i_y_max%nthreads != 0) buffer_division = (image_buffer_size/nthreads) + 1;
    else buffer_division = (image_buffer_size/nthreads);

    for(int i = 0; i < image_buffer_size; i++){
        fwrite(image_collection[i/buffer_division][i%buffer_division], 1 , 3, file);
    };

    fclose(file);
};

void compute_mandelbrot(){

    #pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();

        double z_x;
        double z_y;
        double z_x_squared;
        double z_y_squared;
        double escape_radius_squared = 4;

        int iteration;
        int i_x;
        int i_y;

        double c_x;
        double c_y;

        unsigned char ** image_buffer_thread;

        // Variavel que coordena a divisão do tamanho dos buffers e divisao de trabalho.
        int buffer_division;
        if(i_y_max%nthreads != 0) buffer_division = (image_buffer_size/nthreads) + 1;
        else buffer_division = (image_buffer_size/nthreads);

        // Trecho abaixo realiza allocação nos moldes do antigo allocate_image_buffer()
        // Só que para buffers menores um para cada thread.
        int rgb_size = 3;
        image_buffer_thread = (unsigned char **) malloc(sizeof(unsigned char *) * buffer_division);
        for(int i = 0; i < buffer_division; i++){
            image_buffer_thread[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
        };

        // Ajuste inicial do x, pois uma thread pode comecar onde a ultima parou em um x diferente
        i_x = (buffer_division* tid)%i_y_max;

        for(i_y = tid*(buffer_division/i_y_max); i_y <= ((tid + 1)*buffer_division)/i_y_max && i_y < i_y_max ; i_y++){
            c_y = c_y_min + i_y * pixel_height;

            if(fabs(c_y) < pixel_height / 2) {
                c_y = 0.0;
            };

            // x máximo para termino da divisao de trabalho.
            int x_division;
            if (i_y == ((tid + 1)*buffer_division)/i_y_max) x_division = ((tid + 1)*buffer_division)%i_x_max;
            else x_division = i_x_max ;

            for(; i_x < i_x_max && i_x < x_division; i_x++){
                //printf("Thread %d: pos(%d,%d)\n", tid, i_x, i_y);
                c_x         = c_x_min + i_x * pixel_width;

                z_x         = 0.0;
                z_y         = 0.0;

                z_x_squared = 0.0;
                z_y_squared = 0.0;

                for(iteration = 0;
                    iteration < iteration_max && \
                    ((z_x_squared + z_y_squared) < escape_radius_squared);
                    iteration++){
                    z_y         = 2 * z_x * z_y + c_y;
                    z_x         = z_x_squared - z_y_squared + c_x;

                    z_x_squared = z_x * z_x;
                    z_y_squared = z_y * z_y;
                };

                // Á propria threads realiza o update de buffers
                int color;
                if(iteration == iteration_max){
                    image_buffer_thread[((i_y_max * i_y) + i_x)%buffer_division][0] = colors[gradient_size][0];
                    image_buffer_thread[((i_y_max * i_y) + i_x)%buffer_division][1] = colors[gradient_size][1];
                    image_buffer_thread[((i_y_max * i_y) + i_x)%buffer_division][2] = colors[gradient_size][2];
                }
                else{
                    color = iteration % gradient_size;
                    image_buffer_thread[((i_y_max * i_y) + i_x)%buffer_division][0] = colors[color][0];
                    image_buffer_thread[((i_y_max * i_y) + i_x)%buffer_division][1] = colors[color][1];
                    image_buffer_thread[((i_y_max * i_y) + i_x)%buffer_division][2] = colors[color][2];
                };
            };
            // Somente no inicio de cada thread x pode ser diferente de zero, devido a lugar
            // que a ultima parou.
            i_x = 0;
        };
        //Anexa o buffer da thread no array de buffers.
        image_collection[tid] = image_buffer_thread;
    }
};

int main(int argc, char *argv[]){

    // Recebe numero de threads, caso não há argv[6]: threads = 1
    if(argc > 6) nthreads = atoi(argv[6]);

    init(argc, argv);

    // Cria um array de ponteiros para cada buffer criado ser anexado
    allocate_image_collection();

    compute_mandelbrot();

    write_to_file();

    return 0;
};
