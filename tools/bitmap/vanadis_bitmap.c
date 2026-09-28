/*
 * vanadis_bitmap.c
 *
 * Simple bitmap postprocessor for a structured 2-D Vanadis output plane.
 *
 * Input format (one node per line):
 *     coord1  coord2  concentration
 *
 * The node ordering is expected to be structured and row-major:
 * coord1 varies fastest, coord2 varies between rows.  Grid dimensions are
 * detected automatically from the coordinates, so they do not need to be
 * hard-coded.
 *
 * The renderer uses bilinear Q4 interpolation inside each rectangular cell.
 * Negative interpolated concentrations are clipped to zero for visualization,
 * reproducing the behavior of the original Vanadis bitmap utility.
 *
 * Default invocation reproduces the historical test setup:
 *     vanadis_bitmap
 *
 * Equivalent explicit invocation:
 *     vanadis_bitmap X_02000_001.TXT vanadis_ground.bmp \
 *                    -1500 1500 -1500 1500 500 500
 *
 * Build:
 *     gcc -O2 -std=c11 -Wall -Wextra -pedantic vanadis_bitmap.c \
 *         -o vanadis_bitmap.exe
 */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PALETTE_SIZE 15
#define BMP_HEADER_SIZE 54

#define DEFAULT_INPUT  "X_02000_001.TXT"
#define DEFAULT_OUTPUT "vanadis_ground.bmp"
#define DEFAULT_XMIN   (-1500.0)
#define DEFAULT_XMAX   ( 1500.0)
#define DEFAULT_YMIN   (-1500.0)
#define DEFAULT_YMAX   ( 1500.0)
#define DEFAULT_WIDTH  500
#define DEFAULT_HEIGHT 500

typedef struct {
    double x;
    double y;
    double c;
} Node;

typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
} RGB;

/* Historical Vanadis bitmap palette. */
static const RGB palette[PALETTE_SIZE] = {
    {  0, 200, 215},
    {247, 160,   0},
    {239, 152,   0},
    {231, 144,   0},
    {223, 136,   0},
    {215, 128,   0},
    {207, 120,   0},
    {199, 112,   0},
    {191, 104,   0},
    {175,  88,   0},
    {167,  80,   0},
    {159,  72,   0},
    {151,  64,   0},
    {143,  56,   0},
    {135,  48,   0}
};

static void usage(const char *program)
{
    fprintf(stderr,
        "Usage:\n"
        "  %s\n"
        "  %s INPUT OUTPUT\n"
        "  %s INPUT OUTPUT XMIN XMAX YMIN YMAX WIDTH HEIGHT\n\n"
        "Defaults:\n"
        "  INPUT  = %s\n"
        "  OUTPUT = %s\n"
        "  bounds = [%g, %g] x [%g, %g]\n"
        "  size   = %d x %d pixels\n",
        program, program, program,
        DEFAULT_INPUT, DEFAULT_OUTPUT,
        DEFAULT_XMIN, DEFAULT_XMAX, DEFAULT_YMIN, DEFAULT_YMAX,
        DEFAULT_WIDTH, DEFAULT_HEIGHT);
}

static int nearly_equal(double a, double b)
{
    const double scale = fmax(1.0, fmax(fabs(a), fabs(b)));
    return fabs(a - b) <= 1.0e-10 * scale;
}

static int parse_double(const char *text, double *value)
{
    char *end = NULL;
    errno = 0;
    *value = strtod(text, &end);
    return errno == 0 && end != text && *end == '\0';
}

static int parse_int(const char *text, int *value)
{
    char *end = NULL;
    long v;

    errno = 0;
    v = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || v <= 0 || v > 1000000) {
        return 0;
    }

    *value = (int)v;
    return 1;
}

static int read_nodes(const char *filename, Node **nodes_out, size_t *count_out,
                      double *max_concentration)
{
    FILE *fp = fopen(filename, "r");
    Node *nodes = NULL;
    size_t count = 0;
    size_t capacity = 0;
    double max_c = 0.0;

    if (fp == NULL) {
        fprintf(stderr, "Error: cannot open input file '%s'.\n", filename);
        return 0;
    }

    for (;;) {
        Node node;
        int rc = fscanf(fp, "%lf %lf %lf", &node.x, &node.y, &node.c);

        if (rc == EOF) {
            break;
        }
        if (rc != 3) {
            fprintf(stderr, "Error: malformed input near node %zu in '%s'.\n",
                    count + 1, filename);
            free(nodes);
            fclose(fp);
            return 0;
        }

        if (count == capacity) {
            size_t new_capacity = capacity == 0 ? 1024 : capacity * 2;
            Node *tmp = realloc(nodes, new_capacity * sizeof(*nodes));
            if (tmp == NULL) {
                fprintf(stderr, "Error: out of memory while reading '%s'.\n", filename);
                free(nodes);
                fclose(fp);
                return 0;
            }
            nodes = tmp;
            capacity = new_capacity;
        }

        nodes[count++] = node;
        if (node.c > max_c) {
            max_c = node.c;
        }
    }

    fclose(fp);

    if (count < 4) {
        fprintf(stderr, "Error: '%s' contains too few nodes.\n", filename);
        free(nodes);
        return 0;
    }

    *nodes_out = nodes;
    *count_out = count;
    *max_concentration = max_c;
    return 1;
}

static int detect_grid(const Node *nodes, size_t count,
                       size_t *nx_out, size_t *ny_out,
                       double **x_axis_out, double **y_axis_out)
{
    size_t nx = 1;
    size_t ny;
    double *x_axis = NULL;
    double *y_axis = NULL;

    while (nx < count && nearly_equal(nodes[nx].y, nodes[0].y)) {
        ++nx;
    }

    if (nx < 2 || count % nx != 0) {
        fprintf(stderr,
                "Error: input is not a rectangular structured grid "
                "(could not infer dimensions).\n");
        return 0;
    }

    ny = count / nx;
    if (ny < 2) {
        fprintf(stderr, "Error: detected grid has fewer than two rows.\n");
        return 0;
    }

    x_axis = malloc(nx * sizeof(*x_axis));
    y_axis = malloc(ny * sizeof(*y_axis));
    if (x_axis == NULL || y_axis == NULL) {
        fprintf(stderr, "Error: out of memory while constructing grid axes.\n");
        free(x_axis);
        free(y_axis);
        return 0;
    }

    for (size_t ix = 0; ix < nx; ++ix) {
        x_axis[ix] = nodes[ix].x;
    }
    for (size_t iy = 0; iy < ny; ++iy) {
        y_axis[iy] = nodes[iy * nx].y;
    }

    /* Verify that every row/column uses the same structured coordinates. */
    for (size_t iy = 0; iy < ny; ++iy) {
        for (size_t ix = 0; ix < nx; ++ix) {
            const Node *node = &nodes[iy * nx + ix];
            if (!nearly_equal(node->x, x_axis[ix]) ||
                !nearly_equal(node->y, y_axis[iy])) {
                fprintf(stderr,
                        "Error: node ordering is not a rectangular structured grid "
                        "at row %zu, column %zu.\n",
                        iy + 1, ix + 1);
                free(x_axis);
                free(y_axis);
                return 0;
            }
        }
    }

    *nx_out = nx;
    *ny_out = ny;
    *x_axis_out = x_axis;
    *y_axis_out = y_axis;
    return 1;
}

static int axis_is_monotonic(const double *axis, size_t n)
{
    int direction;

    if (n < 2 || nearly_equal(axis[0], axis[1])) {
        return 0;
    }

    direction = axis[1] > axis[0] ? 1 : -1;
    for (size_t i = 1; i < n; ++i) {
        if (direction > 0) {
            if (!(axis[i] > axis[i - 1])) {
                return 0;
            }
        } else {
            if (!(axis[i] < axis[i - 1])) {
                return 0;
            }
        }
    }
    return 1;
}

/* Return the lower interval index i such that value lies between axis[i]
 * and axis[i+1].  Works for strictly increasing or decreasing axes. */
static int find_interval(const double *axis, size_t n, double value, size_t *index)
{
    const int increasing = axis[n - 1] > axis[0];

    if (increasing) {
        if (value < axis[0] || value > axis[n - 1]) {
            return 0;
        }
        if (nearly_equal(value, axis[n - 1])) {
            *index = n - 2;
            return 1;
        }
        for (size_t i = 0; i + 1 < n; ++i) {
            if (value >= axis[i] && value < axis[i + 1]) {
                *index = i;
                return 1;
            }
        }
    } else {
        if (value > axis[0] || value < axis[n - 1]) {
            return 0;
        }
        if (nearly_equal(value, axis[n - 1])) {
            *index = n - 2;
            return 1;
        }
        for (size_t i = 0; i + 1 < n; ++i) {
            if (value <= axis[i] && value > axis[i + 1]) {
                *index = i;
                return 1;
            }
        }
    }

    return 0;
}

static double interpolate_q4(const Node *nodes, size_t nx,
                             size_t ix, size_t iy,
                             double x, double y)
{
    const Node *n1 = &nodes[iy * nx + ix];
    const Node *n2 = &nodes[iy * nx + (ix + 1)];
    const Node *n3 = &nodes[(iy + 1) * nx + (ix + 1)];
    const Node *n4 = &nodes[(iy + 1) * nx + ix];

    const double dx = n2->x - n1->x;
    const double dy = n4->y - n1->y;

    double ksi;
    double eta;
    double N1, N2, N3, N4;

    if (dx == 0.0 || dy == 0.0) {
        return 0.0;
    }

    ksi = 2.0 * (x - n1->x) / dx - 1.0;
    eta = 2.0 * (y - n1->y) / dy - 1.0;

    N1 = 0.25 * (1.0 - ksi) * (1.0 - eta);
    N2 = 0.25 * (1.0 + ksi) * (1.0 - eta);
    N3 = 0.25 * (1.0 + ksi) * (1.0 + eta);
    N4 = 0.25 * (1.0 - ksi) * (1.0 + eta);

    return N1 * n1->c + N2 * n2->c + N3 * n3->c + N4 * n4->c;
}

static void put_u16_le(uint8_t *p, uint16_t value)
{
    p[0] = (uint8_t)(value & 0xffu);
    p[1] = (uint8_t)((value >> 8) & 0xffu);
}

static void put_u32_le(uint8_t *p, uint32_t value)
{
    p[0] = (uint8_t)(value & 0xffu);
    p[1] = (uint8_t)((value >> 8) & 0xffu);
    p[2] = (uint8_t)((value >> 16) & 0xffu);
    p[3] = (uint8_t)((value >> 24) & 0xffu);
}

static int write_bmp32(const char *filename,
                       const uint8_t *pixels,
                       int width, int height)
{
    uint8_t header[BMP_HEADER_SIZE] = {0};
    const uint64_t image_size_64 = (uint64_t)width * (uint64_t)height * 4u;
    const uint64_t file_size_64 = BMP_HEADER_SIZE + image_size_64;
    FILE *fp;

    if (image_size_64 > UINT32_MAX || file_size_64 > UINT32_MAX) {
        fprintf(stderr, "Error: bitmap is too large for this writer.\n");
        return 0;
    }

    header[0] = 'B';
    header[1] = 'M';
    put_u32_le(&header[2], (uint32_t)file_size_64);
    put_u32_le(&header[10], BMP_HEADER_SIZE);
    put_u32_le(&header[14], 40);                 /* BITMAPINFOHEADER size */
    put_u32_le(&header[18], (uint32_t)width);
    put_u32_le(&header[22], (uint32_t)height);  /* positive = bottom-up */
    put_u16_le(&header[26], 1);                 /* planes */
    put_u16_le(&header[28], 32);                /* BGRA, 32 bits/pixel */
    put_u32_le(&header[34], (uint32_t)image_size_64);

    fp = fopen(filename, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Error: cannot create output file '%s'.\n", filename);
        return 0;
    }

    if (fwrite(header, 1, sizeof(header), fp) != sizeof(header) ||
        fwrite(pixels, 1, (size_t)image_size_64, fp) != (size_t)image_size_64) {
        fprintf(stderr, "Error: failed while writing '%s'.\n", filename);
        fclose(fp);
        return 0;
    }

    if (fclose(fp) != 0) {
        fprintf(stderr, "Error: failed while closing '%s'.\n", filename);
        return 0;
    }

    return 1;
}

int main(int argc, char **argv)
{
    const char *input_file = DEFAULT_INPUT;
    const char *output_file = DEFAULT_OUTPUT;
    double xmin = DEFAULT_XMIN;
    double xmax = DEFAULT_XMAX;
    double ymin = DEFAULT_YMIN;
    double ymax = DEFAULT_YMAX;
    int width = DEFAULT_WIDTH;
    int height = DEFAULT_HEIGHT;

    Node *nodes = NULL;
    size_t node_count = 0;
    size_t nx = 0, ny = 0;
    double *x_axis = NULL;
    double *y_axis = NULL;
    double max_concentration = 0.0;
    uint8_t *pixels = NULL;
    size_t pixel_bytes;

    if (argc == 3) {
        input_file = argv[1];
        output_file = argv[2];
    } else if (argc == 9) {
        input_file = argv[1];
        output_file = argv[2];

        if (!parse_double(argv[3], &xmin) || !parse_double(argv[4], &xmax) ||
            !parse_double(argv[5], &ymin) || !parse_double(argv[6], &ymax) ||
            !parse_int(argv[7], &width) || !parse_int(argv[8], &height)) {
            fprintf(stderr, "Error: invalid numeric command-line argument.\n\n");
            usage(argv[0]);
            return EXIT_FAILURE;
        }
    } else if (argc != 1) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (!(xmin < xmax) || !(ymin < ymax)) {
        fprintf(stderr, "Error: invalid rendering bounds.\n");
        return EXIT_FAILURE;
    }

    if (!read_nodes(input_file, &nodes, &node_count, &max_concentration)) {
        return EXIT_FAILURE;
    }

    if (!detect_grid(nodes, node_count, &nx, &ny, &x_axis, &y_axis)) {
        free(nodes);
        return EXIT_FAILURE;
    }

    if (!axis_is_monotonic(x_axis, nx) || !axis_is_monotonic(y_axis, ny)) {
        fprintf(stderr, "Error: grid axes must be strictly monotonic.\n");
        free(x_axis);
        free(y_axis);
        free(nodes);
        return EXIT_FAILURE;
    }

    if (xmin < fmin(x_axis[0], x_axis[nx - 1]) ||
        xmax > fmax(x_axis[0], x_axis[nx - 1]) ||
        ymin < fmin(y_axis[0], y_axis[ny - 1]) ||
        ymax > fmax(y_axis[0], y_axis[ny - 1])) {
        fprintf(stderr,
                "Error: requested rendering bounds extend outside the input grid.\n"
                "Grid bounds: x=[%.10g, %.10g], y=[%.10g, %.10g]\n",
                fmin(x_axis[0], x_axis[nx - 1]),
                fmax(x_axis[0], x_axis[nx - 1]),
                fmin(y_axis[0], y_axis[ny - 1]),
                fmax(y_axis[0], y_axis[ny - 1]));
        free(x_axis);
        free(y_axis);
        free(nodes);
        return EXIT_FAILURE;
    }

    if (max_concentration <= 0.0) {
        fprintf(stderr,
                "Warning: maximum nodal concentration is <= 0; "
                "bitmap will contain only the background color.\n");
    }

    if ((uint64_t)width * (uint64_t)height > SIZE_MAX / 4u) {
        fprintf(stderr, "Error: requested bitmap dimensions are too large.\n");
        free(x_axis);
        free(y_axis);
        free(nodes);
        return EXIT_FAILURE;
    }

    pixel_bytes = (size_t)width * (size_t)height * 4u;
    pixels = malloc(pixel_bytes);
    if (pixels == NULL) {
        fprintf(stderr, "Error: cannot allocate %zu bytes for bitmap pixels.\n",
                pixel_bytes);
        free(x_axis);
        free(y_axis);
        free(nodes);
        return EXIT_FAILURE;
    }

    for (int row = 0; row < height; ++row) {
        /* Preserve the historical mapping: divide by height/width rather than
         * height-1/width-1.  With a positive BMP height, row 0 is the bottom. */
        const double y = ymin + (ymax - ymin) * (double)row / (double)height;
        size_t iy;

        if (!find_interval(y_axis, ny, y, &iy)) {
            fprintf(stderr, "Error: cannot locate y=%.10g in grid.\n", y);
            free(pixels);
            free(x_axis);
            free(y_axis);
            free(nodes);
            return EXIT_FAILURE;
        }

        for (int col = 0; col < width; ++col) {
            const double x = xmin + (xmax - xmin) * (double)col / (double)width;
            size_t ix;
            double concentration;
            int color_index = 0;
            size_t p;

            if (!find_interval(x_axis, nx, x, &ix)) {
                fprintf(stderr, "Error: cannot locate x=%.10g in grid.\n", x);
                free(pixels);
                free(x_axis);
                free(y_axis);
                free(nodes);
                return EXIT_FAILURE;
            }

            concentration = interpolate_q4(nodes, nx, ix, iy, x, y);

            /* Visualization choice retained from the original program. */
            if (concentration < 0.0) {
                concentration = 0.0;
            }

            if (max_concentration > 0.0) {
                const double step = max_concentration / (double)(PALETTE_SIZE - 1);
                color_index = (int)floor(concentration / step);
                if (color_index < 0) {
                    color_index = 0;
                }
                if (color_index >= PALETTE_SIZE) {
                    color_index = PALETTE_SIZE - 1;
                }
            }

            p = ((size_t)row * (size_t)width + (size_t)col) * 4u;
            pixels[p + 0] = palette[color_index].b;
            pixels[p + 1] = palette[color_index].g;
            pixels[p + 2] = palette[color_index].r;
            pixels[p + 3] = 255;
        }
    }

    if (!write_bmp32(output_file, pixels, width, height)) {
        free(pixels);
        free(x_axis);
        free(y_axis);
        free(nodes);
        return EXIT_FAILURE;
    }

    printf("Input:              %s\n", input_file);
    printf("Detected grid:      nx=%zu, ny=%zu (%zu nodes)\n", nx, ny, node_count);
    printf("Maximum nodal C:    %.10g\n", max_concentration);
    printf("Rendering bounds:   x=[%.10g, %.10g], y=[%.10g, %.10g]\n",
           xmin, xmax, ymin, ymax);
    printf("Bitmap size:        %d x %d\n", width, height);
    printf("Output:             %s\n", output_file);

    free(pixels);
    free(x_axis);
    free(y_axis);
    free(nodes);
    return EXIT_SUCCESS;
}
