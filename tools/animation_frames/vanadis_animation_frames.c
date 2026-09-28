/*
   bmp_xyz_corrected.c

   Jeden program do generowania BMP dla:
       X_?????_???.TXT
       Y_?????_???.TXT
       Z_?????_???.TXT

   Geometria i kolejnosc danych wynikaja BEZPOSREDNIO z kodu Fortran:

     nax = 7
     nay = 20
     naz = 20

   siatka:
     do jz = 1,naz
        do jy = 1,nay
           do jx = 1,nax
              xx = dx(jx)
              yy = dy(jy)
              zz = dz(jz)

   --------------------------------------------------------------------------
   X - stala wspolrzedna xx(jx):
       plik zawiera: YY, ZZ, S
       20 x 20 wezlow
       YY zmienia sie najszybciej.

   Y - stala wspolrzedna yy(jy):
       plik zawiera: ZZ, XX, S
       20 x 7 wezlow
       XX (7 wysokosci) zmienia sie najszybciej.

   Z - stala wspolrzedna zz(jz):
       plik zawiera: YY, XX, S
       20 x 7 wezlow
       XX (7 wysokosci) zmienia sie najszybciej.

   --------------------------------------------------------------------------
   Skala:

   GLOBAL:
       - jedno wspolne maksimum dla wszystkich plikow X, Y, Z
         i wszystkich czasow.

   FAMILY:
       - jedna wspolna skala dla tej samej litery plaszczyzny
         i tego samego indeksu przez wszystkie czasy.
       - np. X_*_001.TXT to jedna family,
             X_*_002.TXT to druga family,
             Y_*_001.TXT to trzecia family.

   LOCAL:
       - osobne maksimum i osobna skala dla kazdego pliku.

   Uzycie:
       bmp_xyz.exe C:\tmp\output
       bmp_xyz.exe C:\tmp\output global
       bmp_xyz.exe C:\tmp\output family
       bmp_xyz.exe C:\tmp\output local

   Dodatkowo:
       - wartosci wspolrzednych sa wypisywane przy osiach,
       - piksel < 0.1% maksimum danego pliku jest ustawiany na zero,
       - na osiach symetrycznych pokazujemy min, 0, max,
       - dla przekrojow pionowych: poziomo -600,0,600,
         pionowo 0...420,
       - zachowana jest jednakowa skala metr/piksel na obu osiach.

   Kompilacja:
       gcc bmp_xyz_axes.c -O3 -o bmp_xyz.exe -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <ctype.h>

#ifdef _WIN32
    #include <windows.h>
    #define PATH_SEP '\\'
#else
    #include <dirent.h>
    #define PATH_SEP '/'
#endif

#define COLORS 15

#define NAX 7
#define NAY 20
#define NAZ 20

/*
   X ma zachowywac sie jak dzialajacy bmp_with_global_scale.c.
*/
#define X_PLOT_WIDTH   500
#define X_PLOT_HEIGHT  500

#define X_VIEW_MIN    -1500.0
#define X_VIEW_MAX     1500.0

/*
   Przekroje Y/Z:
       poziomo = Y lub Z
       pionowo = X = wysokosc 0..420

   Wycinamy centralny zakres -600..600 m z pelnej domeny -4000..4000.
   Zachowujemy fizyczna proporcje osi:
       1200 m -> 500 px
       420 m  -> 175 px
*/
#define V_VIEW_MIN       -600.0
#define V_VIEW_MAX        600.0
#define V_HEIGHT_MIN        0.0
#define V_HEIGHT_MAX      420.0

#define V_PLOT_WIDTH      500
#define V_PIXEL_PER_METER ((double)V_PLOT_WIDTH / (V_VIEW_MAX - V_VIEW_MIN))
#define V_PLOT_HEIGHT     ((int)((V_HEIGHT_MAX - V_HEIGHT_MIN) * V_PIXEL_PER_METER + 0.5))

#define MARGIN_LEFT   100
#define MARGIN_TOP     25
#define MARGIN_RIGHT   20
#define MARGIN_BOTTOM  50

#define LEGEND_GAP     30
#define LEGEND_BAR_W   30
#define LEGEND_TICK_W   8
#define LABEL_SCALE     2
#define LABEL_GAP       8
#define LABEL_AREA_W  120

typedef enum {
    SCALE_GLOBAL,
    SCALE_FAMILY,
    SCALE_LOCAL
} ScaleMode;

typedef struct {
    char name[512];
    char family;                /* X, Y albo Z */
    unsigned long long time_sec;
    unsigned int plane_index;
    double local_max;
    double family_max;
} InputFile;

typedef struct {
    int width;
    int height;
    unsigned char *pixels;      /* BGRA */
} Image;


/* ------------------------------------------------------------------------- */
/* Osie dokladnie z wsp_na_o() */

static const double DX[NAX] = {
      0.0,
     20.0,
     35.0,
     85.0,
    135.0,
    230.0,
    420.0
};

static const double DY[NAY] = {
    -4000.0,
    -1500.0,
     -900.0,
     -700.0,
     -500.0,
     -400.0,
     -270.0,
     -180.0,
      -90.0,
      -25.0,
       25.0,
       90.0,
      180.0,
      270.0,
      400.0,
      500.0,
      700.0,
      900.0,
     1500.0,
     4000.0
};

static const double DZ[NAZ] = {
    -4000.0,
    -1500.0,
     -900.0,
     -700.0,
     -500.0,
     -400.0,
     -270.0,
     -180.0,
      -90.0,
      -25.0,
       25.0,
       90.0,
      180.0,
      270.0,
      400.0,
      500.0,
      700.0,
      900.0,
     1500.0,
     4000.0
};


/* ------------------------------------------------------------------------- */
/* Prototypy */

static int parse_input_name(const char *name, InputFile *out);
static int cmp_input_files(const void *a, const void *b);

static void build_path(char *dst, size_t dst_size,
                       const char *dir, const char *name);

static int collect_files(const char *dir,
                         InputFile **files_out,
                         size_t *count_out);

static int scan_file_max(const char *dir, InputFile *f);

static int locate_interval(const double *axis, int n, double value);

static int process_file(const char *dir,
                        const InputFile *f,
                        double scale_max);

static int render_X(const char *input_path,
                    Image *img,
                    double scale_max,
                    double pixel_cutoff,
                    const int *R,
                    const int *G,
                    const int *B);

static int render_vertical(const char *input_path,
                           Image *img,
                           char family,
                           double scale_max,
                           double pixel_cutoff,
                           const int *R,
                           const int *G,
                           const int *B);

static void init_palette(int *R, int *G, int *B);
static int concentration_to_color(double s, double scale_max);

static int image_init(Image *img,
                      int width,
                      int height,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a);

static void image_free(Image *img);

static void set_pixel(Image *img,
                      int x,
                      int y,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a);

static void fill_rect(Image *img,
                      int x,
                      int y,
                      int w,
                      int h,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a);

static void draw_rect(Image *img,
                      int x,
                      int y,
                      int w,
                      int h,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a);

static const unsigned char *get_glyph(char c);

static void draw_char(Image *img,
                      int x,
                      int y,
                      char c,
                      int scale,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a);

static void draw_text(Image *img,
                      int x,
                      int y,
                      const char *text,
                      int scale,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a);

static void format_value(char *buf, size_t buf_size, double value);

static int text_width_px(const char *text, int scale);

static void draw_axis_dimensions(Image *img,
                                 int plot_x0,
                                 int plot_y0,
                                 int plot_width,
                                 int plot_height,
                                 double x_min,
                                 double x_max,
                                 double y_min,
                                 double y_max);

static void draw_legend(Image *img,
                        int plot_x0,
                        int plot_y0,
                        int plot_width,
                        int plot_height,
                        double scale_max,
                        const int *R,
                        const int *G,
                        const int *B);

static int write_bmp32_topdown(const char *path, const Image *img);


/* ------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
    const char *dir = (argc >= 2) ? argv[1] : ".";

    ScaleMode scale_mode = SCALE_GLOBAL;

    if (argc >= 3) {
        if (strcmp(argv[2], "global") == 0) {
            scale_mode = SCALE_GLOBAL;
        }
        else if (strcmp(argv[2], "family") == 0) {
            scale_mode = SCALE_FAMILY;
        }
        else if (strcmp(argv[2], "local") == 0) {
            scale_mode = SCALE_LOCAL;
        }
        else {
            fprintf(stderr,
                    "Nieznany tryb skali: %s\n"
                    "Dozwolone: global, family, local\n",
                    argv[2]);
            return 1;
        }
    }

    InputFile *files = NULL;
    size_t count = 0;

    if (collect_files(dir, &files, &count) != 0) {
        fprintf(stderr,
                "Blad podczas odczytu katalogu: %s\n",
                dir);
        return 1;
    }

    if (count == 0) {
        printf("Nie znaleziono plikow X/Y/Z_?????_???.TXT "
               "w katalogu: %s\n",
               dir);
        free(files);
        return 0;
    }

    qsort(files, count, sizeof(InputFile), cmp_input_files);

    /*
       Etap 1.
       Obliczamy maksimum dodatnie z kazdego pliku.
    */

    double max_global = 0.0;

    printf("Znaleziono %zu plikow X/Y/Z.\n", count);
    printf("Wyznaczanie maksimow stezenia...\n");

    for (size_t i = 0; i < count; ++i) {

        if (scan_file_max(dir, &files[i]) != 0) {
            free(files);
            return 1;
        }

        if (files[i].local_max > max_global)
            max_global = files[i].local_max;
    }

    if (max_global <= 0.0)
        max_global = 1.0;

    /*
       FAMILY:
       ta sama litera plaszczyzny + ten sam indeks plaszczyzny,
       przez wszystkie czasy.

       Przyklady osobnych rodzin:
           X_*_001.TXT
           X_*_002.TXT
           Y_*_001.TXT
           Z_*_001.TXT
    */
    for (size_t i = 0; i < count; ++i) {

        double family_max = 0.0;

        for (size_t j = 0; j < count; ++j) {

            if (files[j].family == files[i].family &&
                files[j].plane_index == files[i].plane_index &&
                files[j].local_max > family_max) {

                family_max = files[j].local_max;
            }
        }

        if (family_max <= 0.0)
            family_max = 1.0;

        files[i].family_max = family_max;
    }

    printf("\nMaksimum global = %.12g\n", max_global);

    if (scale_mode == SCALE_LOCAL)
        printf("Tryb skali: LOCAL = osobna skala dla kazdego pliku\n");
    else if (scale_mode == SCALE_FAMILY)
        printf("Tryb skali: FAMILY = wspolna skala dla tej samej litery i indeksu plaszczyzny przez wszystkie czasy\n");
    else
        printf("Tryb skali: GLOBAL = jedna skala dla wszystkich X/Y/Z i wszystkich czasow\n");

    printf("\n");

    /*
       Etap 2.
       Generowanie BMP.
    */

    size_t ok = 0;

    for (size_t i = 0; i < count; ++i) {

        double scale_max = 1.0;

        if (scale_mode == SCALE_LOCAL) {

            scale_max = files[i].local_max;

            if (scale_max <= 0.0)
                scale_max = 1.0;
        }
        else if (scale_mode == SCALE_FAMILY) {

            scale_max = files[i].family_max;
        }
        else {

            scale_max = max_global;
        }

        printf("[%zu/%zu] %c  t=%llu s  plane=%03u  "
               "local_max=%.6g  scale=%.6g\n",
               i + 1,
               count,
               files[i].family,
               files[i].time_sec,
               files[i].plane_index,
               files[i].local_max,
               scale_max);

        printf("          prog 0.1%% local_max = %.6g\n",
               0.001 * files[i].local_max);

        if (process_file(dir,
                         &files[i],
                         scale_max) == 0) {
            ++ok;
        }
        else {
            fprintf(stderr,
                    "  -> BLAD dla %s\n",
                    files[i].name);
        }
    }

    printf("\nGotowe: %zu/%zu plikow BMP.\n",
           ok,
           count);

    free(files);

    return (ok == count) ? 0 : 2;
}


/* ------------------------------------------------------------------------- */
/*
   Nazwy:
       X_00100_001.TXT
       Y_00100_010.TXT
       Z_00100_010.TXT
*/

static int parse_input_name(const char *name, InputFile *out)
{
    if (strlen(name) != 15)
        return 0;

    char family =
        (char)toupper((unsigned char)name[0]);

    if (family != 'X' &&
        family != 'Y' &&
        family != 'Z')
        return 0;

    if (name[1] != '_' ||
        name[7] != '_' ||
        name[11] != '.')
        return 0;

    for (int i = 2; i <= 6; ++i) {
        if (!isdigit((unsigned char)name[i]))
            return 0;
    }

    for (int i = 8; i <= 10; ++i) {
        if (!isdigit((unsigned char)name[i]))
            return 0;
    }

    if (tolower((unsigned char)name[12]) != 't' ||
        tolower((unsigned char)name[13]) != 'x' ||
        tolower((unsigned char)name[14]) != 't')
        return 0;

    char time_buf[6];
    char plane_buf[4];

    memcpy(time_buf, name + 2, 5);
    time_buf[5] = '\0';

    memcpy(plane_buf, name + 8, 3);
    plane_buf[3] = '\0';

    memset(out, 0, sizeof(*out));

    snprintf(out->name,
             sizeof(out->name),
             "%s",
             name);

    out->family = family;

    out->time_sec =
        strtoull(time_buf, NULL, 10);

    out->plane_index =
        (unsigned int)strtoul(plane_buf, NULL, 10);

    out->local_max = 0.0;
    out->family_max = 0.0;

    return 1;
}


/* ------------------------------------------------------------------------- */

static int cmp_input_files(const void *a, const void *b)
{
    const InputFile *fa =
        (const InputFile *)a;

    const InputFile *fb =
        (const InputFile *)b;

    if (fa->family < fb->family) return -1;
    if (fa->family > fb->family) return 1;

    if (fa->time_sec < fb->time_sec) return -1;
    if (fa->time_sec > fb->time_sec) return 1;

    if (fa->plane_index < fb->plane_index) return -1;
    if (fa->plane_index > fb->plane_index) return 1;

    return strcmp(fa->name, fb->name);
}


/* ------------------------------------------------------------------------- */

static void build_path(char *dst,
                       size_t dst_size,
                       const char *dir,
                       const char *name)
{
    size_t n = strlen(dir);

    if (n > 0 &&
        (dir[n - 1] == '/' ||
         dir[n - 1] == '\\')) {

        snprintf(dst,
                 dst_size,
                 "%s%s",
                 dir,
                 name);
    }
    else {

        snprintf(dst,
                 dst_size,
                 "%s%c%s",
                 dir,
                 PATH_SEP,
                 name);
    }
}


/* ------------------------------------------------------------------------- */

static int add_input_file(InputFile **files,
                          size_t *count,
                          size_t *capacity,
                          const InputFile *f)
{
    if (*count == *capacity) {

        size_t new_capacity =
            (*capacity == 0)
            ? 64
            : (*capacity * 2);

        InputFile *tmp =
            (InputFile *)realloc(
                *files,
                new_capacity * sizeof(InputFile));

        if (!tmp)
            return -1;

        *files = tmp;
        *capacity = new_capacity;
    }

    (*files)[*count] = *f;

    ++(*count);

    return 0;
}


/* ------------------------------------------------------------------------- */

static int collect_files(const char *dir,
                         InputFile **files_out,
                         size_t *count_out)
{
    InputFile *files = NULL;

    size_t count = 0;
    size_t capacity = 0;

#ifdef _WIN32

    char pattern[1024];

    size_t n = strlen(dir);

    if (n > 0 &&
        (dir[n - 1] == '/' ||
         dir[n - 1] == '\\')) {

        snprintf(pattern,
                 sizeof(pattern),
                 "%s*_*_*.*",
                 dir);
    }
    else {

        snprintf(pattern,
                 sizeof(pattern),
                 "%s\\*_*_*.*",
                 dir);
    }

    WIN32_FIND_DATAA fd;

    HANDLE hFind =
        FindFirstFileA(pattern, &fd);

    if (hFind == INVALID_HANDLE_VALUE) {

        DWORD err = GetLastError();

        if (err == ERROR_FILE_NOT_FOUND) {

            *files_out = NULL;
            *count_out = 0;

            return 0;
        }

        return -1;
    }

    do {

        if (!(fd.dwFileAttributes &
              FILE_ATTRIBUTE_DIRECTORY)) {

            InputFile f;

            if (parse_input_name(fd.cFileName,
                                 &f)) {

                if (add_input_file(&files,
                                   &count,
                                   &capacity,
                                   &f) != 0) {

                    FindClose(hFind);
                    free(files);
                    return -1;
                }
            }
        }

    } while (FindNextFileA(hFind, &fd));

    FindClose(hFind);

#else

    DIR *dp = opendir(dir);

    if (!dp)
        return -1;

    struct dirent *entry;

    while ((entry = readdir(dp)) != NULL) {

        InputFile f;

        if (parse_input_name(entry->d_name,
                             &f)) {

            if (add_input_file(&files,
                               &count,
                               &capacity,
                               &f) != 0) {

                closedir(dp);
                free(files);
                return -1;
            }
        }
    }

    closedir(dp);

#endif

    *files_out = files;
    *count_out = count;

    return 0;
}


/* ------------------------------------------------------------------------- */
/*
   Skanuje tylko trzecia kolumne.

   X:
       20*20 = 400 rekordow

   Y/Z:
       20*7 = 140 rekordow
*/

static int scan_file_max(const char *dir, InputFile *f)
{
    char input_path[1024];

    build_path(input_path,
               sizeof(input_path),
               dir,
               f->name);

    FILE *fp =
        fopen(input_path, "r");

    if (!fp) {
        perror(input_path);
        return -1;
    }

    int node_count;

    if (f->family == 'X')
        node_count = NAY * NAZ;
    else
        node_count = NAX * NAZ;

    double max_value = 0.0;

    for (int i = 0;
         i < node_count;
         ++i) {

        double a, b, s;

        if (fscanf(fp,
                   "%lf %lf %lf",
                   &a,
                   &b,
                   &s) != 3) {

            fprintf(stderr,
                    "Niepoprawny format %s, "
                    "rekord %d/%d\n",
                    input_path,
                    i + 1,
                    node_count);

            fclose(fp);
            return -1;
        }

        if (s > max_value)
            max_value = s;
    }

    fclose(fp);

    f->local_max = max_value;

    return 0;
}


/* ------------------------------------------------------------------------- */
/* O( log N ) */

static int locate_interval(const double *axis,
                           int n,
                           double value)
{
    if (n < 2)
        return -1;

    if (value < axis[0] ||
        value > axis[n - 1])
        return -1;

    if (value == axis[n - 1])
        return n - 2;

    int lo = 0;
    int hi = n - 2;

    while (lo <= hi) {

        int mid =
            lo + (hi - lo) / 2;

        if (value < axis[mid]) {

            hi = mid - 1;
        }
        else if (value >= axis[mid + 1]) {

            lo = mid + 1;
        }
        else {

            return mid;
        }
    }

    return -1;
}


/* ------------------------------------------------------------------------- */

static void init_palette(int *R, int *G, int *B)
{
    R[0] = 0;   G[0] = 200; B[0] = 215;

    R[1] = 247; G[1] = 160; B[1] = 0;
    R[2] = 239; G[2] = 152; B[2] = 0;
    R[3] = 231; G[3] = 144; B[3] = 0;
    R[4] = 223; G[4] = 136; B[4] = 0;
    R[5] = 215; G[5] = 128; B[5] = 0;
    R[6] = 207; G[6] = 120; B[6] = 0;
    R[7] = 199; G[7] = 112; B[7] = 0;
    R[8] = 191; G[8] = 104; B[8] = 0;
    R[9] = 175; G[9] = 88;  B[9] = 0;

    R[10] = 167; G[10] = 80; B[10] = 0;
    R[11] = 159; G[11] = 72; B[11] = 0;
    R[12] = 151; G[12] = 64; B[12] = 0;
    R[13] = 143; G[13] = 56; B[13] = 0;
    R[14] = 135; G[14] = 48; B[14] = 0;
}


/* ------------------------------------------------------------------------- */
/*
   Zachowanie jak w dzialajacym programie X:

       ss < 0 -> 0
       kolor = floor(ss / step)

   indeks 0 = turkusowe tlo
*/

static int concentration_to_color(double s,
                                  double scale_max)
{
    if (s <= 0.0)
        return 0;

    double step =
        scale_max / (COLORS - 1);

    if (step <= 0.0)
        return 0;

    int idx =
        (int)fabs(s / step);

    if (idx < 0)
        idx = 0;

    if (idx >= COLORS)
        idx = COLORS - 1;

    return idx;
}


/* ------------------------------------------------------------------------- */
/* Obraz */

static int image_init(Image *img,
                      int width,
                      int height,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a)
{
    img->width = width;
    img->height = height;

    img->pixels =
        (unsigned char *)malloc(
            (size_t)width *
            height *
            4);

    if (!img->pixels)
        return -1;

    for (int y = 0;
         y < height;
         ++y) {

        for (int x = 0;
             x < width;
             ++x) {

            size_t p =
                ((size_t)y *
                 width +
                 x) *
                4;

            img->pixels[p + 0] = b;
            img->pixels[p + 1] = g;
            img->pixels[p + 2] = r;
            img->pixels[p + 3] = a;
        }
    }

    return 0;
}


/* ------------------------------------------------------------------------- */

static void image_free(Image *img)
{
    free(img->pixels);

    img->pixels = NULL;
    img->width = 0;
    img->height = 0;
}


/* ------------------------------------------------------------------------- */

static void set_pixel(Image *img,
                      int x,
                      int y,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a)
{
    if (x < 0 ||
        x >= img->width ||
        y < 0 ||
        y >= img->height)
        return;

    size_t p =
        ((size_t)y *
         img->width +
         x) *
        4;

    img->pixels[p + 0] = b;
    img->pixels[p + 1] = g;
    img->pixels[p + 2] = r;
    img->pixels[p + 3] = a;
}


/* ------------------------------------------------------------------------- */

static void fill_rect(Image *img,
                      int x,
                      int y,
                      int w,
                      int h,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a)
{
    if (w <= 0 || h <= 0)
        return;

    int x0 =
        (x < 0) ? 0 : x;

    int y0 =
        (y < 0) ? 0 : y;

    int x1 =
        (x + w > img->width)
        ? img->width
        : x + w;

    int y1 =
        (y + h > img->height)
        ? img->height
        : y + h;

    for (int yy = y0;
         yy < y1;
         ++yy) {

        for (int xx = x0;
             xx < x1;
             ++xx) {

            set_pixel(img,
                      xx,
                      yy,
                      r, g, b, a);
        }
    }
}


/* ------------------------------------------------------------------------- */

static void draw_rect(Image *img,
                      int x,
                      int y,
                      int w,
                      int h,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a)
{
    for (int xx = x;
         xx < x + w;
         ++xx) {

        set_pixel(img,
                  xx,
                  y,
                  r, g, b, a);

        set_pixel(img,
                  xx,
                  y + h - 1,
                  r, g, b, a);
    }

    for (int yy = y;
         yy < y + h;
         ++yy) {

        set_pixel(img,
                  x,
                  yy,
                  r, g, b, a);

        set_pixel(img,
                  x + w - 1,
                  yy,
                  r, g, b, a);
    }
}


/* ------------------------------------------------------------------------- */
/* Czcionka 5x7 */

static const unsigned char GLYPH_0[7] = {14,17,19,21,25,17,14};
static const unsigned char GLYPH_1[7] = {4,12,4,4,4,4,14};
static const unsigned char GLYPH_2[7] = {14,17,1,2,4,8,31};
static const unsigned char GLYPH_3[7] = {30,1,1,14,1,1,30};
static const unsigned char GLYPH_4[7] = {2,6,10,18,31,2,2};
static const unsigned char GLYPH_5[7] = {31,16,16,30,1,1,30};
static const unsigned char GLYPH_6[7] = {14,16,16,30,17,17,14};
static const unsigned char GLYPH_7[7] = {31,1,2,4,8,8,8};
static const unsigned char GLYPH_8[7] = {14,17,17,14,17,17,14};
static const unsigned char GLYPH_9[7] = {14,17,17,15,1,1,14};

static const unsigned char GLYPH_DOT[7] =
    {0,0,0,0,0,12,12};

static const unsigned char GLYPH_MINUS[7] =
    {0,0,0,31,0,0,0};

static const unsigned char GLYPH_PLUS[7] =
    {0,4,4,31,4,4,0};

static const unsigned char GLYPH_E[7] =
    {31,16,16,30,16,16,31};

static const unsigned char GLYPH_e[7] =
    {0,0,14,17,31,16,14};

static const unsigned char GLYPH_SPACE[7] =
    {0,0,0,0,0,0,0};


/* ------------------------------------------------------------------------- */

static const unsigned char *get_glyph(char c)
{
    switch (c) {

        case '0': return GLYPH_0;
        case '1': return GLYPH_1;
        case '2': return GLYPH_2;
        case '3': return GLYPH_3;
        case '4': return GLYPH_4;
        case '5': return GLYPH_5;
        case '6': return GLYPH_6;
        case '7': return GLYPH_7;
        case '8': return GLYPH_8;
        case '9': return GLYPH_9;

        case '.': return GLYPH_DOT;
        case '-': return GLYPH_MINUS;
        case '+': return GLYPH_PLUS;

        case 'E': return GLYPH_E;
        case 'e': return GLYPH_e;

        default:
            return GLYPH_SPACE;
    }
}


/* ------------------------------------------------------------------------- */

static void draw_char(Image *img,
                      int x,
                      int y,
                      char c,
                      int scale,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a)
{
    const unsigned char *glyph =
        get_glyph(c);

    for (int row = 0;
         row < 7;
         ++row) {

        for (int col = 0;
             col < 5;
             ++col) {

            if (glyph[row] &
                (1 << (4 - col))) {

                fill_rect(img,
                          x + col * scale,
                          y + row * scale,
                          scale,
                          scale,
                          r, g, b, a);
            }
        }
    }
}


/* ------------------------------------------------------------------------- */

static void draw_text(Image *img,
                      int x,
                      int y,
                      const char *text,
                      int scale,
                      unsigned char r,
                      unsigned char g,
                      unsigned char b,
                      unsigned char a)
{
    int cursor = x;

    while (*text) {

        draw_char(img,
                  cursor,
                  y,
                  *text,
                  scale,
                  r, g, b, a);

        cursor +=
            6 * scale;

        ++text;
    }
}


/* ------------------------------------------------------------------------- */

static void format_value(char *buf,
                         size_t buf_size,
                         double value)
{
    snprintf(buf,
             buf_size,
             "%.3g",
             value);
}


/* ------------------------------------------------------------------------- */
/* Szerokosc tekstu dla czcionki 5x7: 5 pikseli + 1 piksel odstepu. */

static int text_width_px(const char *text, int scale)
{
    if (!text || !*text)
        return 0;

    int n = (int)strlen(text);

    return n * 6 * scale - scale;
}


/* ------------------------------------------------------------------------- */
/*
   Opisy wspolrzednych na osiach.

   X:
       - wartosc minimalna na lewym koncu,
       - 0 w miejscu wynikajacym z geometrii, jesli 0 jest w zakresie,
       - wartosc maksymalna na prawym koncu.

   Y:
       - minimum na dole,
       - 0 jesli lezy wewnatrz zakresu,
       - maksimum na gorze.

   Rysujemy rowniez male kreski (ticki).
*/

static void draw_axis_dimensions(Image *img,
                                 int plot_x0,
                                 int plot_y0,
                                 int plot_width,
                                 int plot_height,
                                 double x_min,
                                 double x_max,
                                 double y_min,
                                 double y_max)
{
    const int scale = LABEL_SCALE;
    const int tick = 6;
    const int text_gap = 5;

    char buf[64];

    /* ------------------------------------------------------------- */
    /* OS POZIOMA */

    /* min */
    format_value(buf, sizeof(buf), x_min);

    draw_text(img,
              plot_x0,
              plot_y0 + plot_height + tick + text_gap,
              buf,
              scale,
              0,0,0,255);

    for (int t = 0; t < tick; ++t) {
        set_pixel(img,
                  plot_x0,
                  plot_y0 + plot_height + t,
                  0,0,0,255);
    }

    /* max */
    format_value(buf, sizeof(buf), x_max);

    int tw = text_width_px(buf, scale);

    draw_text(img,
              plot_x0 + plot_width - tw,
              plot_y0 + plot_height + tick + text_gap,
              buf,
              scale,
              0,0,0,255);

    for (int t = 0; t < tick; ++t) {
        set_pixel(img,
                  plot_x0 + plot_width - 1,
                  plot_y0 + plot_height + t,
                  0,0,0,255);
    }

    /* zero, jesli lezy wewnatrz zakresu */
    if (x_min < 0.0 && x_max > 0.0) {

        double frac =
            (0.0 - x_min) /
            (x_max - x_min);

        int px =
            plot_x0 +
            (int)lround(
                frac *
                (plot_width - 1));

        strcpy(buf, "0");

        tw = text_width_px(buf, scale);

        draw_text(img,
                  px - tw / 2,
                  plot_y0 + plot_height + tick + text_gap,
                  buf,
                  scale,
                  0,0,0,255);

        for (int t = 0; t < tick; ++t) {
            set_pixel(img,
                      px,
                      plot_y0 + plot_height + t,
                      0,0,0,255);
        }
    }


    /* ------------------------------------------------------------- */
    /* OS PIONOWA */

    /* minimum - dol */
    format_value(buf, sizeof(buf), y_min);

    tw = text_width_px(buf, scale);

    draw_text(img,
              plot_x0 - tick - text_gap - tw,
              plot_y0 + plot_height - 7 * scale,
              buf,
              scale,
              0,0,0,255);

    for (int t = 0; t < tick; ++t) {
        set_pixel(img,
                  plot_x0 - t - 1,
                  plot_y0 + plot_height - 1,
                  0,0,0,255);
    }

    /* maksimum - gora */
    format_value(buf, sizeof(buf), y_max);

    tw = text_width_px(buf, scale);

    draw_text(img,
              plot_x0 - tick - text_gap - tw,
              plot_y0,
              buf,
              scale,
              0,0,0,255);

    for (int t = 0; t < tick; ++t) {
        set_pixel(img,
                  plot_x0 - t - 1,
                  plot_y0,
                  0,0,0,255);
    }

    /*
       0 na osi pionowej tylko wtedy, gdy jest W SRODKU zakresu.
       Gdy y_min == 0, opis minimum juz pokazuje zero.
    */
    if (y_min < 0.0 && y_max > 0.0) {

        double frac =
            (0.0 - y_min) /
            (y_max - y_min);

        /*
           W obrazie Y rosnie w dol, wiec odwracamy.
        */
        int py =
            plot_y0 +
            plot_height -
            1 -
            (int)lround(
                frac *
                (plot_height - 1));

        strcpy(buf, "0");

        tw = text_width_px(buf, scale);

        draw_text(img,
                  plot_x0 - tick - text_gap - tw,
                  py - (7 * scale) / 2,
                  buf,
                  scale,
                  0,0,0,255);

        for (int t = 0; t < tick; ++t) {
            set_pixel(img,
                      plot_x0 - t - 1,
                      py,
                      0,0,0,255);
        }
    }
}


/* ------------------------------------------------------------------------- */

static void draw_legend(Image *img,
                        int plot_x0,
                        int plot_y0,
                        int plot_width,
                        int plot_height,
                        double scale_max,
                        const int *R,
                        const int *G,
                        const int *B)
{
    int legend_x =
        plot_x0 +
        plot_width +
        LEGEND_GAP;

    int legend_y =
        plot_y0;

    for (int y = 0;
         y < plot_height;
         ++y) {

        double t =
            1.0 -
            (double)y /
            (double)(plot_height - 1);

        int idx =
            (int)floor(
                t *
                (COLORS - 1) +
                0.5);

        if (idx < 0)
            idx = 0;

        if (idx >= COLORS)
            idx = COLORS - 1;

        for (int x = 0;
             x < LEGEND_BAR_W;
             ++x) {

            set_pixel(img,
                      legend_x + x,
                      legend_y + y,
                      (unsigned char)R[idx],
                      (unsigned char)G[idx],
                      (unsigned char)B[idx],
                      255);
        }
    }

    draw_rect(img,
              legend_x,
              legend_y,
              LEGEND_BAR_W,
              plot_height,
              0,0,0,255);

    /*
       7 opisow:
       0, 1/6, 2/6, ... max
    */

    for (int k = 0;
         k <= 6;
         ++k) {

        double frac =
            (double)k / 6.0;

        double value =
            frac * scale_max;

        int y =
            legend_y +
            plot_height -
            1 -
            (int)lround(
                frac *
                (plot_height - 1));

        for (int t = 0;
             t < LEGEND_TICK_W;
             ++t) {

            set_pixel(img,
                      legend_x +
                      LEGEND_BAR_W +
                      t,
                      y,
                      0,0,0,255);
        }

        char buf[64];

        format_value(buf,
                     sizeof(buf),
                     value);

        draw_text(img,
                  legend_x +
                  LEGEND_BAR_W +
                  LEGEND_TICK_W +
                  LABEL_GAP,
                  y -
                  (7 * LABEL_SCALE) / 2,
                  buf,
                  LABEL_SCALE,
                  0,0,0,255);
    }
}


/* ------------------------------------------------------------------------- */
/*
   X - KOD ZACHOWUJACY LOGIKE DZIALAJACEGO
       bmp_with_global_scale.c

   Dane:
       col1 = Y
       col2 = Z
       col3 = S

   Kolejnosc:
       dla jednego Z:
           Y(1)..Y(20)
*/

static int render_X(const char *input_path,
                    Image *img,
                    double scale_max,
                    double pixel_cutoff,
                    const int *R,
                    const int *G,
                    const int *B)
{
    const int node_count =
        NAY * NAZ;

    double c1[node_count + 1];
    double c2[node_count + 1];
    double s [node_count + 1];

    FILE *fp =
        fopen(input_path, "r");

    if (!fp) {
        perror(input_path);
        return -1;
    }

    for (int i = 1;
         i <= node_count;
         ++i) {

        if (fscanf(fp,
                   "%lf %lf %lf",
                   &c1[i],
                   &c2[i],
                   &s[i]) != 3) {

            fprintf(stderr,
                    "Niepoprawny format %s, "
                    "rekord %d/%d\n",
                    input_path,
                    i,
                    node_count);

            fclose(fp);
            return -1;
        }
    }

    fclose(fp);

    /*
       Os 1 = pierwsze 20 rekordow.
       Os 2 = pierwszy rekord kazdego bloku 20.
    */

    double axis1[NAY];
    double axis2[NAZ];

    for (int i = 0;
         i < NAY;
         ++i) {

        axis1[i] =
            c1[i + 1];
    }

    for (int i = 0;
         i < NAZ;
         ++i) {

        axis2[i] =
            c2[i * NAY + 1];
    }

    /*
       Prekalkulacja jak w dzialajacym programie:
       zakres -1500...1500.
    */

    int e1_for_x[X_PLOT_WIDTH];
    double ksi_for_x[X_PLOT_WIDTH];

    for (int ii = 0;
         ii < X_PLOT_WIDTH;
         ++ii) {

        double x =
            X_VIEW_MIN +
            (X_VIEW_MAX - X_VIEW_MIN) /
            (double)X_PLOT_WIDTH *
            (double)ii;

        int e1 =
            locate_interval(axis1,
                            NAY,
                            x);

        e1_for_x[ii] =
            e1;

        ksi_for_x[ii] =
            0.0;

        if (e1 >= 0) {

            double left =
                axis1[e1];

            double right =
                axis1[e1 + 1];

            double d =
                right - left;

            if (fabs(d) > 1e-15) {

                ksi_for_x[ii] =
                    2.0 *
                    (x - left) /
                    d -
                    1.0;
            }
        }
    }

    int e2_for_y[X_PLOT_HEIGHT];
    double eta_for_y[X_PLOT_HEIGHT];

    for (int jj = 0;
         jj < X_PLOT_HEIGHT;
         ++jj) {

        double y =
            X_VIEW_MIN +
            (X_VIEW_MAX - X_VIEW_MIN) /
            (double)X_PLOT_HEIGHT *
            (double)(
                X_PLOT_HEIGHT -
                1 -
                jj);

        int e2 =
            locate_interval(axis2,
                            NAZ,
                            y);

        e2_for_y[jj] =
            e2;

        eta_for_y[jj] =
            0.0;

        if (e2 >= 0) {

            double bottom =
                axis2[e2];

            double top =
                axis2[e2 + 1];

            double d =
                top - bottom;

            if (fabs(d) > 1e-15) {

                eta_for_y[jj] =
                    2.0 *
                    (y - bottom) /
                    d -
                    1.0;
            }
        }
    }

    const int plot_x0 =
        MARGIN_LEFT;

    const int plot_y0 =
        MARGIN_TOP;

    for (int jj = 0;
         jj < X_PLOT_HEIGHT;
         ++jj) {

        int e2 =
            e2_for_y[jj];

        double eta =
            eta_for_y[jj];

        double eta_minus =
            1.0 - eta;

        double eta_plus =
            1.0 + eta;

        for (int ii = 0;
             ii < X_PLOT_WIDTH;
             ++ii) {

            int e1 =
                e1_for_x[ii];

            int RGB_PIX =
                0;

            if (e1 >= 0 &&
                e2 >= 0) {

                double ksi =
                    ksi_for_x[ii];

                double ksi_minus =
                    1.0 - ksi;

                double ksi_plus =
                    1.0 + ksi;

                /*
                   Dokladnie taki sam uklad jak
                   w bmp_with_global_scale:
                */

                int n1 =
                    e2 * NAY +
                    e1 +
                    1;

                int n2 =
                    n1 +
                    1;

                int n4 =
                    n1 +
                    NAY;

                int n3 =
                    n4 +
                    1;

                double f1 =
                    0.25 *
                    ksi_minus *
                    eta_minus;

                double f2 =
                    0.25 *
                    ksi_plus *
                    eta_minus;

                double f3 =
                    0.25 *
                    ksi_plus *
                    eta_plus;

                double f4 =
                    0.25 *
                    ksi_minus *
                    eta_plus;

                double ss =
                    f1 * s[n1] +
                    f2 * s[n2] +
                    f3 * s[n3] +
                    f4 * s[n4];

                if (ss < 0.0)
                    ss = 0.0;

                /*
                   Prog wizualizacji:
                   jezeli interpolowana wartosc piksela jest mniejsza
                   niz 0.1% maksimum TEGO KONKRETNEGO PLIKU,
                   traktujemy ja jako zero.

                   Prog nie zalezy od trybu skali kolorow
                   (global/family/local).
                */
                if (ss < pixel_cutoff)
                    ss = 0.0;

                RGB_PIX =
                    concentration_to_color(
                        ss,
                        scale_max);
            }

            set_pixel(img,
                      plot_x0 + ii,
                      plot_y0 + jj,
                      (unsigned char)R[RGB_PIX],
                      (unsigned char)G[RGB_PIX],
                      (unsigned char)B[RGB_PIX],
                      255);
        }
    }

    draw_rect(img,
              plot_x0,
              plot_y0,
              X_PLOT_WIDTH,
              X_PLOT_HEIGHT,
              0,0,0,255);

    draw_axis_dimensions(img,
                         plot_x0,
                         plot_y0,
                         X_PLOT_WIDTH,
                         X_PLOT_HEIGHT,
                         X_VIEW_MIN,
                         X_VIEW_MAX,
                         X_VIEW_MIN,
                         X_VIEW_MAX);

    draw_legend(img,
                plot_x0,
                plot_y0,
                X_PLOT_WIDTH,
                X_PLOT_HEIGHT,
                scale_max,
                R,G,B);

    return 0;
}


/* ------------------------------------------------------------------------- */
/*
   Y oraz Z.

   Z FORTRANA:

   Y:
       write x3(...), x1(...), s(...)
       czyli:
           col1 = Z
           col2 = X
           col3 = S

   Z:
       write x2(...), x1(...), s(...)
       czyli:
           col1 = Y
           col2 = X
           col3 = S

   W obu przypadkach:
       20 pozycji poziomych
       x 7 wysokosci.

   Kolejnosc:
       horizontal(1), X(1..7)
       horizontal(2), X(1..7)
       ...
       horizontal(20), X(1..7)
*/

static int render_vertical(const char *input_path,
                           Image *img,
                           char family,
                           double scale_max,
                           double pixel_cutoff,
                           const int *R,
                           const int *G,
                           const int *B)
{
    const int node_count =
        NAX * NAY;

    /*
       grid[horizontal][height]
    */

    double grid[NAY][NAX];

    double horizontal_from_file[NAY];
    double height_from_file[NAX];

    FILE *fp =
        fopen(input_path, "r");

    if (!fp) {
        perror(input_path);
        return -1;
    }

    for (int ih = 0;
         ih < NAY;
         ++ih) {

        for (int ix = 0;
             ix < NAX;
             ++ix) {

            double horizontal;
            double height;
            double s;

            if (fscanf(fp,
                       "%lf %lf %lf",
                       &horizontal,
                       &height,
                       &s) != 3) {

                fprintf(stderr,
                        "Niepoprawny format %s, "
                        "rekord %d/%d\n",
                        input_path,
                        ih * NAX +
                        ix +
                        1,
                        node_count);

                fclose(fp);
                return -1;
            }

            /*
               Stężenie fizycznie nie powinno być ujemne.
               Ujemne wartości w danych są oscylacją numeryczną.
               Dla Y/Z obcinamy je W WĘZŁACH przed interpolacją,
               żeby wartości ujemne nie wycinały sztucznych dziur
               w chmurze stężenia.
            */
            if (s < 0.0)
                s = 0.0;

            grid[ih][ix] =
                s;

            if (ix == 0)
                horizontal_from_file[ih] =
                    horizontal;

            if (ih == 0)
                height_from_file[ix] =
                    height;
        }
    }

    fclose(fp);

    /*
       Kontrola kolejnosci danych.

       Y -> pierwsza kolumna powinna zgadzac sie z DZ.
       Z -> pierwsza kolumna powinna zgadzac sie z DY.

       W praktyce DY == DZ, ale kontrola pozostaje jawna.
    */

    const double *horizontal_axis =
        (family == 'Y')
        ? DZ
        : DY;

    for (int i = 0;
         i < NAY;
         ++i) {

        if (fabs(
            horizontal_from_file[i] -
            horizontal_axis[i]) >
            1e-6) {

            fprintf(stderr,
                    "UWAGA: %s, horizontal[%d]=%.12g, "
                    "oczekiwano %.12g\n",
                    input_path,
                    i,
                    horizontal_from_file[i],
                    horizontal_axis[i]);
        }
    }

    for (int i = 0;
         i < NAX;
         ++i) {

        if (fabs(
            height_from_file[i] -
            DX[i]) >
            1e-6) {

            fprintf(stderr,
                    "UWAGA: %s, height[%d]=%.12g, "
                    "oczekiwano %.12g\n",
                    input_path,
                    i,
                    height_from_file[i],
                    DX[i]);
        }
    }

    /*
       Prekalkulacja osi poziomej.
    */

    int eh_for_x[V_PLOT_WIDTH];
    double ksi_for_x[V_PLOT_WIDTH];

    for (int px = 0;
         px < V_PLOT_WIDTH;
         ++px) {

        double horizontal =
            V_VIEW_MIN +
            (V_VIEW_MAX - V_VIEW_MIN) /
            (double)V_PLOT_WIDTH *
            (double)px;

        int eh =
            locate_interval(
                horizontal_axis,
                NAY,
                horizontal);

        eh_for_x[px] =
            eh;

        ksi_for_x[px] =
            0.0;

        if (eh >= 0) {

            double left =
                horizontal_axis[eh];

            double right =
                horizontal_axis[eh + 1];

            double d =
                right - left;

            if (fabs(d) > 1e-15) {

                ksi_for_x[px] =
                    2.0 *
                    (horizontal -
                     left) /
                    d -
                    1.0;
            }
        }
    }

    /*
       Prekalkulacja wysokosci.

       dol = X=0
       gora = X=420
    */

    int ex_for_y[V_PLOT_HEIGHT];
    double eta_for_y[V_PLOT_HEIGHT];

    for (int py = 0;
         py < V_PLOT_HEIGHT;
         ++py) {

        double height =
            V_HEIGHT_MIN +
            (V_HEIGHT_MAX -
             V_HEIGHT_MIN) /
            (double)V_PLOT_HEIGHT *
            (double)(
                V_PLOT_HEIGHT -
                1 -
                py);

        int ex =
            locate_interval(
                DX,
                NAX,
                height);

        ex_for_y[py] =
            ex;

        eta_for_y[py] =
            0.0;

        if (ex >= 0) {

            double bottom =
                DX[ex];

            double top =
                DX[ex + 1];

            double d =
                top - bottom;

            if (fabs(d) > 1e-15) {

                eta_for_y[py] =
                    2.0 *
                    (height -
                     bottom) /
                    d -
                    1.0;
            }
        }
    }

    const int plot_x0 =
        MARGIN_LEFT;

    const int plot_y0 =
        MARGIN_TOP;

    /*
       Interpolacja biliniowa.

       grid[horizontal][height]

          n4 -------- n3
           |           |
           |           |
          n1 -------- n2

       tutaj:
          n1 = grid[eh    ][ex    ]
          n2 = grid[eh + 1][ex    ]
          n3 = grid[eh + 1][ex + 1]
          n4 = grid[eh    ][ex + 1]
    */

    for (int py = 0;
         py < V_PLOT_HEIGHT;
         ++py) {

        int ex =
            ex_for_y[py];

        double eta =
            eta_for_y[py];

        double eta_minus =
            1.0 - eta;

        double eta_plus =
            1.0 + eta;

        for (int px = 0;
             px < V_PLOT_WIDTH;
             ++px) {

            int eh =
                eh_for_x[px];

            int RGB_PIX =
                0;

            if (eh >= 0 &&
                ex >= 0) {

                double ksi =
                    ksi_for_x[px];

                double ksi_minus =
                    1.0 - ksi;

                double ksi_plus =
                    1.0 + ksi;

                double f1 =
                    0.25 *
                    ksi_minus *
                    eta_minus;

                double f2 =
                    0.25 *
                    ksi_plus *
                    eta_minus;

                double f3 =
                    0.25 *
                    ksi_plus *
                    eta_plus;

                double f4 =
                    0.25 *
                    ksi_minus *
                    eta_plus;

                double ss =
                    f1 *
                    grid[eh][ex] +

                    f2 *
                    grid[eh + 1][ex] +

                    f3 *
                    grid[eh + 1][ex + 1] +

                    f4 *
                    grid[eh][ex + 1];

                if (ss < 0.0)
                    ss = 0.0;

                /*
                   Prog wizualizacji:
                   jezeli interpolowana wartosc piksela jest mniejsza
                   niz 0.1% maksimum TEGO KONKRETNEGO PLIKU,
                   traktujemy ja jako zero.

                   Prog nie zalezy od trybu skali kolorow
                   (global/family/local).
                */
                if (ss < pixel_cutoff)
                    ss = 0.0;

                RGB_PIX =
                    concentration_to_color(
                        ss,
                        scale_max);
            }

            set_pixel(img,
                      plot_x0 + px,
                      plot_y0 + py,
                      (unsigned char)R[RGB_PIX],
                      (unsigned char)G[RGB_PIX],
                      (unsigned char)B[RGB_PIX],
                      255);
        }
    }

    draw_rect(img,
              plot_x0,
              plot_y0,
              V_PLOT_WIDTH,
              V_PLOT_HEIGHT,
              0,0,0,255);

    draw_axis_dimensions(img,
                         plot_x0,
                         plot_y0,
                         V_PLOT_WIDTH,
                         V_PLOT_HEIGHT,
                         V_VIEW_MIN,
                         V_VIEW_MAX,
                         V_HEIGHT_MIN,
                         V_HEIGHT_MAX);

    draw_legend(img,
                plot_x0,
                plot_y0,
                V_PLOT_WIDTH,
                V_PLOT_HEIGHT,
                scale_max,
                R,G,B);

    return 0;
}


/* ------------------------------------------------------------------------- */

static int process_file(const char *dir,
                        const InputFile *f,
                        double scale_max)
{
    int R[COLORS];
    int G[COLORS];
    int B[COLORS];

    init_palette(R,G,B);

    int plot_width;
    int plot_height;

    if (f->family == 'X') {

        plot_width =
            X_PLOT_WIDTH;

        plot_height =
            X_PLOT_HEIGHT;
    }
    else {

        plot_width =
            V_PLOT_WIDTH;

        plot_height =
            V_PLOT_HEIGHT;
    }

    int image_width =
        MARGIN_LEFT +
        plot_width +
        LEGEND_GAP +
        LEGEND_BAR_W +
        LEGEND_TICK_W +
        LABEL_GAP +
        LABEL_AREA_W +
        MARGIN_RIGHT;

    int image_height =
        MARGIN_TOP +
        plot_height +
        MARGIN_BOTTOM;

    Image img;

    if (image_init(&img,
                   image_width,
                   image_height,
                   255,255,255,255) != 0) {

        fprintf(stderr,
                "Brak pamieci dla %s\n",
                f->name);

        return -1;
    }

    char input_path[1024];

    build_path(input_path,
               sizeof(input_path),
               dir,
               f->name);

    int result;

    /*
       0.1% maksimum tego konkretnego obrazu/pliku.
       Przyklad:
           local_max = 12.0  -> cutoff = 0.012
    */
    double pixel_cutoff =
        0.001 * f->local_max;

    if (pixel_cutoff < 0.0)
        pixel_cutoff = 0.0;

    if (f->family == 'X') {

        result =
            render_X(input_path,
                     &img,
                     scale_max,
                     pixel_cutoff,
                     R,G,B);
    }
    else {

        result =
            render_vertical(input_path,
                            &img,
                            f->family,
                            scale_max,
                            pixel_cutoff,
                            R,G,B);
    }

    if (result != 0) {

        image_free(&img);
        return -1;
    }

    char output_name[512];

    snprintf(output_name,
             sizeof(output_name),
             "%s",
             f->name);

    char *dot =
        strrchr(output_name, '.');

    if (dot)
        strcpy(dot, ".bmp");

    char output_path[1024];

    build_path(output_path,
               sizeof(output_path),
               dir,
               output_name);

    if (write_bmp32_topdown(
            output_path,
            &img) != 0) {

        image_free(&img);
        return -1;
    }

    image_free(&img);

    printf("  -> %s\n",
           output_name);

    return 0;
}


/* ------------------------------------------------------------------------- */
/* BMP 32-bit top-down */

static int write_bmp32_topdown(const char *path,
                               const Image *img)
{
    FILE *fout =
        fopen(path, "wb");

    if (!fout) {
        perror(path);
        return -1;
    }

    const uint32_t pixel_bytes =
        (uint32_t)(
            (size_t)img->width *
            img->height *
            4);

    const uint32_t file_size =
        54u +
        pixel_bytes;

    unsigned char header[54] =
        {0};

    header[0] = 'B';
    header[1] = 'M';

    memcpy(&header[2],
           &file_size,
           4);

    {
        uint32_t pixel_offset = 54;
        uint32_t dib_size = 40;

        int32_t w =
            img->width;

        /*
           Ujemna wysokosc =
           obraz zapisany top-down.
        */
        int32_t h =
            -img->height;

        uint16_t planes = 1;
        uint16_t bpp = 32;

        uint32_t compression = 0;

        uint32_t image_size =
            pixel_bytes;

        memcpy(&header[10],
               &pixel_offset,
               4);

        memcpy(&header[14],
               &dib_size,
               4);

        memcpy(&header[18],
               &w,
               4);

        memcpy(&header[22],
               &h,
               4);

        memcpy(&header[26],
               &planes,
               2);

        memcpy(&header[28],
               &bpp,
               2);

        memcpy(&header[30],
               &compression,
               4);

        memcpy(&header[34],
               &image_size,
               4);
    }

    if (fwrite(header,
               1,
               sizeof(header),
               fout) != sizeof(header) ||

        fwrite(img->pixels,
               1,
               pixel_bytes,
               fout) != pixel_bytes) {

        fprintf(stderr,
                "Blad zapisu pliku %s\n",
                path);

        fclose(fout);

        return -1;
    }

    fclose(fout);

    return 0;
}
