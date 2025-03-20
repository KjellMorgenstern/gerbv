/*
 * gEDA - GNU Electronic Design Automation
 * This file is a part of gerbv.
 *
 *   Copyright (C) 2000-2003 Stefan Petersen (spe@stacken.kth.se)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111 USA
 */

/** \file pick-and-place.c
    \brief PNP (pick-and-place) parsing functions
    \ingroup libgerbv
*/

#include "gerbv.h"

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "gerber.h"
#include "common.h"
#include "csv.h"
#include "pick-and-place.h"

/**
 * Removes quotation marks from a string.
 * 
 * @param str The string to process (will be modified in-place)
 * @return Pointer to the original string
 */
static char* 
pnp_remove_quotes(char* str) {
    if (!str) return str;
    
    char* p = str;
    char* q = str;
    
    while (*p) {
        if (*p != '"') {
            *q++ = *p;
        }
        p++;
    }
    *q = '\0';
    
    return str;
}

static gerbv_net_t* pnp_new_net(gerbv_net_t* net);
static void         pnp_reset_bbox(gerbv_net_t* net);
static void         pnp_init_net(
            gerbv_net_t* net, gerbv_image_t* image, const char* label, gerbv_aperture_state_t apert_state,
            gerbv_interpolation_t interpol
        );
static int          custom_parse_comma_header(
            char* input, char* result_buffer, char* fields[], int max_fields
        );
static gboolean     pnp_parse_header_line(
            char* buf, char* buf0, int* designator_col, int* footprint_col, int* mid_x_col, int* mid_y_col,
            int* ref_x_col, int* ref_y_col, int* pad_x_col, int* pad_y_col, int* layer_col, int* rotation_col,
            int* comment_col, char* delimiter
        );

void
gerb_transf_free(gerbv_transf_t* transf) {
    g_free(transf);
}

void
gerb_transf_reset(gerbv_transf_t* transf) {
    memset(transf, 0, sizeof(gerbv_transf_t));

    transf->r_mat[0][0] = transf->r_mat[1][1] = 1.0; /*off-diagonals 0 diagonals 1 */
    // transf->r_mat[1][0] = transf->r_mat[0][1] = 0.0;
    transf->scale = 1.0;
    // transf->offset[0] = transf->offset[1] = 0.0;

} /* gerb_transf_reset */

gerbv_transf_t*
gerb_transf_new(void) {
    gerbv_transf_t* transf;

    transf = g_new(gerbv_transf_t, 1);
    gerb_transf_reset(transf);

    return transf;
} /* gerb_transf_new */

//! Rotation
/*! append rotation to transformation.
@param transf transformation to be modified
@param angle in rad (counterclockwise rotation) */

void
gerb_transf_rotate(gerbv_transf_t* transf, double angle) {
    double m[2][2];
    double s = sin(angle), c = cos(angle);

    memcpy(m, transf->r_mat, sizeof(m));
    transf->r_mat[0][0] = c * m[0][0] - s * m[1][0];
    transf->r_mat[0][1] = c * m[0][1] - s * m[1][1];
    transf->r_mat[1][0] = s * m[0][0] + c * m[1][0];
    transf->r_mat[1][1] = s * m[0][1] + c * m[1][1];
    //    transf->offset[0] = transf->offset[1] = 0.0; CHECK ME

} /* gerb_transf_rotate */

//! Translation
/*! append translation to transformation.
@param transf transformation to be modified
@param shift_x translation in x direction
@param shift_y translation in y direction */

void
gerb_transf_shift(gerbv_transf_t* transf, double shift_x, double shift_y) {

    transf->offset[0] += shift_x;
    transf->offset[1] += shift_y;

} /* gerb_transf_shift */

void
gerb_transf_apply(double x, double y, gerbv_transf_t* transf, double* out_x, double* out_y) {

    //    x += transf->offset[0];
    //    y += transf->offset[1];
    *out_x = (x * transf->r_mat[0][0] + y * transf->r_mat[0][1]) * transf->scale;
    *out_y = (x * transf->r_mat[1][0] + y * transf->r_mat[1][1]) * transf->scale;
    *out_x += transf->offset[0];
    *out_y += transf->offset[1];

} /* gerb_transf_apply */

void
pick_and_place_reset_bounding_box(gerbv_net_t* net) {
    net->boundingBox.left   = -HUGE_VAL;
    net->boundingBox.right  = HUGE_VAL;
    net->boundingBox.bottom = -HUGE_VAL;
    net->boundingBox.top    = HUGE_VAL;
}

/* Parses a string representing float number with a unit.
 * Default unit can be specified with def_unit. */
static double
pick_and_place_get_float_unit(const char* str, const char* def_unit) {
    double x            = 0.0;
    char   unit_str[41] = {
          0,
    };
    const char* unit = unit_str;

    /* float, optional space, optional unit mm,cm,in,mil */
    sscanf(str, "%lf %40s", &x, unit_str);

    if (unit_str[0] == '\0')
        unit = def_unit;

    /* NOTE: in order of comparability,
     * i.e. "mm" before "m", as "m" will match "mm" */
    if (strstr(unit, "mm")) {
        x /= 25.4;
    } else if (strstr(unit, "in")) {
        /* NOTE: "in" is without scaling. */
    } else if (strstr(unit, "cmil")) {
        x /= 1e5;
    } else if (strstr(unit, "dmil")) {
        x /= 1e4;
    } else if (strstr(unit, "mil")) {
        x /= 1e3;
    } else if (strstr(unit, "km")) {
        x /= 25.4 / 1e6;
    } else if (strstr(unit, "dm")) {
        x /= 25.4 / 100;
    } else if (strstr(unit, "cm")) {
        x /= 25.4 / 10;
    } else if (strstr(unit, "um")) {
        x /= 25.4 * 1e3;
    } else if (strstr(unit, "nm")) {
        x /= 25.4 * 1e6;
    } else if (strstr(unit, "m")) {
        x /= 25.4 / 1e3;
    } else { /* default to "mil" */
        x /= 1e3;
    }

    return x;
} /* pick_and_place_get_float_unit */

/** search a string for a delimiter.
 Must occur at least n times. */
int
pick_and_place_screen_for_delimiter(char* str, int n) {
    char* ptr;
    char  delimiter[4] = "|,;:";
    int   counter[4];
    int   idx, idx_max = 0;

    memset(counter, 0, sizeof(counter));
    for (ptr = str; *ptr; ptr++) {
        switch (*ptr) {
            case '|': idx = 0; break;
            case ',': idx = 1; break;
            case ';': idx = 2; break;
            case ':': idx = 3; break;
            default: continue; break;
        }
        counter[idx]++;
        if (counter[idx] > counter[idx_max]) {
            idx_max = idx;
        }
    }

    if (counter[idx_max] > n) {
        return (unsigned char)delimiter[idx_max];
    } else {
        return -1;
    }
} /* pick_and_place_screen_for_delimiter */

/**Parses the PNP data.
   two lists are filled with the row data.\n One for the scrollable list in the search and select parts interface, the
   other one a mere two columned list, which drives the autocompletion when entering a search.\n It also tries to
   determine the shape of a part and sets  pnp_state->shape accordingly which will be used when drawing the selections
   as an overlay on screen.
   @return the initial node of the pnp_state netlist
 */

GArray*
pick_and_place_parse_file(gerb_file_t* fd) {
    PnpPartData pnpPartData;
    memset(&pnpPartData, 0, sizeof(PnpPartData));
    int   lineCounter = 0, parsedLines = 0;
    int   ret;
    char* row[12];
    char  buf[MAXL + 2], buf0[MAXL + 2];
    char  def_unit[41] = {0};
    double          tmp_x = 0.05, tmp_y = 0.05; /* Default to reasonable values (~100mil total) */
    gerbv_transf_t* tr_rot            = gerb_transf_new();
    GArray*         pnpParseDataArray = g_array_new(FALSE, FALSE, sizeof(PnpPartData));
    gboolean        foundValidDataRow = FALSE;
    const char*     def_unit_prefix = "# X,Y in ";

    /* Column indices - default values if header parsing fails */
    int designator_col = 0;
    int footprint_col = 1;
    int mid_x_col = 2;
    int mid_y_col = 3;
    int ref_x_col = 4;
    int ref_y_col = 5;
    int pad_x_col = 6;
    int pad_y_col = 7;
    int layer_col = 8;
    int rotation_col = 9;
    int comment_col = 10;

    /*
     * many locales redefine "." as "," and so on, so sscanf has problems when
     * reading Pick and Place files using %f format
     */
    setlocale(LC_NUMERIC, "C");

    /* Default to comma as delimiter, but we'll detect it */
    char delimiter = '\0';

    /* Clear buf and buf0 to avoid any garbage data */
    memset(buf, 0, MAXL + 2);
    memset(buf0, 0, MAXL + 2);

    /* Main processing loop */
    while (fgets(buf, MAXL, fd->fd) != NULL) {
        int len = strlen(buf) - 1;
        int i_length = 0, i_width = 0;

        lineCounter += 1; /* next line */

        /* Trim trailing newlines */
        if (len >= 0 && buf[len] == '\n') buf[len--] = 0;
        if (len >= 0 && buf[len] == '\r') buf[len--] = 0;

        /* Check for unit declaration lines */
        if (0 == strncmp(buf, def_unit_prefix, strlen(def_unit_prefix))) {
            sscanf(&buf[strlen(def_unit_prefix)], "%40s.", def_unit);
            printf("Debug: Detected unit declaration: %s\n", def_unit);
            continue;
        } else if (strstr(buf, "X,Y in mils") || strstr(buf, "X, Y in mils")) {
            strcpy(def_unit, "mil");
            printf("Debug: Detected 'mils' as the unit\n");
            continue;
        } else if (strstr(buf, "X,Y in mm") || strstr(buf, "X, Y in mm")) {
            strcpy(def_unit, "mm");
            printf("Debug: Detected 'mm' as the unit\n");
            continue;
        }

        /* Skip empty or short lines */
        if (len <= 5) continue;

        /* Skip special commands that indicate this is not a PnP file */
        if ((buf[0] == '%') || 
            (strncmp(buf, "G54 ", 4) == 0) || 
            (strncmp(buf, "G04 ", 4) == 0)) {
            g_array_free(pnpParseDataArray, TRUE);
            return NULL;
        }

        /* Try to detect and parse header lines */
        if (pnp_check_and_parse_header(buf, lineCounter, buf0, &designator_col, &footprint_col, 
                                      &mid_x_col, &mid_y_col, &ref_x_col, &ref_y_col, 
                                      &pad_x_col, &pad_y_col, &layer_col, &rotation_col, 
                                      &comment_col, &delimiter, &foundValidDataRow)) {
            /* Successfully parsed a header - continue to next line */
            continue;
        }

        /* Skip comment lines */
        if ((buf[0] == '#' || buf[0] == '*')) {
            /* Count these as valid lines so they don't affect our valid percentage */
            parsedLines += 1;
            continue;
        }

        if(!foundValidDataRow) {
            continue;
        }

        /* Parse data row */
        memset(row, 0, sizeof(char*) * 12);
        memset(buf0, 0, MAXL + 2);

        /* Use appropriate parser based on delimiter */
        if (delimiter == ',') {
            ret = custom_parse_comma_header(buf, buf0, row, 11);
        } else {
            ret = csv_row_parse(buf, MAXL, buf0, MAXL, row, 11, (int)delimiter, CSV_QUOTES);
            if (ret <= 0) {
                ret = custom_parse_comma_header(buf, buf0, row, 11);
            }
        }

        if (ret <= 0) {
            printf("Debug: Failed to parse row: %s\n", buf);
            continue;
        }

        /* Skip description rows that might follow headers */
        if (strstr(buf, "Description:") || strstr(buf, "File generated")) {
            printf("Debug: Skipping description/info row\n");
            continue;
        }

        /* Debug output to show what we parsed */
        for (int i = 0; i < ret && i < 11; i++) {
            if (row[i] != NULL) {
                printf("Debug: Data[%d] = '%s'\n", i, row[i]);
            }
        }

        /* Process data based on identified column structure */
        if (ret >= 3) {
            /* Parse designator */
            if (designator_col < ret && row[designator_col] != NULL) {
                snprintf(pnpPartData.designator, sizeof(pnpPartData.designator) - 1, "%s", row[designator_col]);
            } else if (ret > 0 && row[0] != NULL) {
                /* Fallback to first column if designator column is invalid */
                snprintf(pnpPartData.designator, sizeof(pnpPartData.designator) - 1, "%s", row[0]);
            } else {
                continue; /* Skip rows without a designator */
            }

            /* Parse footprint if available */
            if (footprint_col < ret && row[footprint_col] != NULL) {
                snprintf(pnpPartData.footprint, sizeof(pnpPartData.footprint) - 1, "%s", row[footprint_col]);
            } else {
                pnpPartData.footprint[0] = '\0';
            }

            /* Parse layer */
            if (layer_col < ret && row[layer_col] != NULL) {
                snprintf(pnpPartData.layer, sizeof(pnpPartData.layer) - 1, "%s", row[layer_col]);
            } else {
                strcpy(pnpPartData.layer, "TOP"); /* Default to TOP */
            }

            /* Parse comment if available */
            if (comment_col < ret && row[comment_col] != NULL) {
                if (!g_utf8_validate(row[comment_col], -1, NULL)) {
                    gchar* str = g_convert(row[comment_col], strlen(row[comment_col]), 
                                          "UTF-8", "ISO-8859-1", NULL, NULL, NULL);
                    if (str != NULL) {
                        snprintf(pnpPartData.comment, sizeof(pnpPartData.comment) - 1, "%s", str);
                        g_free(str);
                    } else {
                        pnpPartData.comment[0] = '\0';
                    }
                } else {
                    snprintf(pnpPartData.comment, sizeof(pnpPartData.comment) - 1, "%s", row[comment_col]);
                }
            } else {
                pnpPartData.comment[0] = '\0';
            }

            /* Parse X coordinate */
            if (mid_x_col < ret && row[mid_x_col] != NULL) {
                /* Clean up the X coordinate - handle quoted values */
                char x_str[50] = {0};
                strncpy(x_str, row[mid_x_col], sizeof(x_str)-1);

                /* Remove any quotes */
                pnp_remove_quotes(x_str);

                /* Parse X coordinate */
                pnpPartData.mid_x = pick_and_place_get_float_unit(x_str, def_unit);
            } else {
                pnpPartData.mid_x = 0.0;
            }

            /* Parse Y coordinate */
            if (mid_y_col < ret && row[mid_y_col] != NULL) {
                /* Clean up the Y coordinate - handle quoted values */
                char y_str[50] = {0};
                strncpy(y_str, row[mid_y_col], sizeof(y_str)-1);

                /* Remove any quotes */
                pnp_remove_quotes(y_str);

                /* Parse Y coordinate */
                pnpPartData.mid_y = pick_and_place_get_float_unit(y_str, def_unit);
            } else {
                pnpPartData.mid_y = 0.0;
            }

            /* Skip if coordinates are at origin (likely a header row) */
            if ((fabs(pnpPartData.mid_x) < 0.001) && (fabs(pnpPartData.mid_y) < 0.001)) {
                continue;
            }

            /* Parse reference points if available */
            if (ref_x_col < ret && row[ref_x_col] != NULL) {
                pnpPartData.ref_x = pick_and_place_get_float_unit(row[ref_x_col], def_unit);
            } else {
                pnpPartData.ref_x = pnpPartData.mid_x;
            }

            if (ref_y_col < ret && row[ref_y_col] != NULL) {
                pnpPartData.ref_y = pick_and_place_get_float_unit(row[ref_y_col], def_unit);
            } else {
                pnpPartData.ref_y = pnpPartData.mid_y;
            }

            /* Parse pad positions if available */
            if (pad_x_col < ret && row[pad_x_col] != NULL) {
                pnpPartData.pad_x = pick_and_place_get_float_unit(row[pad_x_col], def_unit);
            } else {
                pnpPartData.pad_x = pnpPartData.mid_x + 0.03;
            }

            if (pad_y_col < ret && row[pad_y_col] != NULL) {
                pnpPartData.pad_y = pick_and_place_get_float_unit(row[pad_y_col], def_unit);
            } else {
                pnpPartData.pad_y = pnpPartData.mid_y + 0.03;
            }

            /* Normalize layer name */
            if (strcasecmp(pnpPartData.layer, "top") == 0 || 
                strcasecmp(pnpPartData.layer, "t") == 0 ||
                strcasecmp(pnpPartData.layer, "1") == 0) {
                strcpy(pnpPartData.layer, "TOP");
            } else if (strcasecmp(pnpPartData.layer, "bottom") == 0 || 
                       strcasecmp(pnpPartData.layer, "bot") == 0 || 
                       strcasecmp(pnpPartData.layer, "b") == 0 ||
                       strcasecmp(pnpPartData.layer, "2") == 0) {
                strcpy(pnpPartData.layer, "BOTTOM");
            } else if (strcasecmp(pnpPartData.layer, "top/bottom") == 0 ||
                       strcasecmp(pnpPartData.layer, "both") == 0) {
                strcpy(pnpPartData.layer, "TOP");
            }

            /* Parse rotation */
            pnpPartData.rotation = 0.0;
            if (rotation_col < ret && row[rotation_col] != NULL) {
                char rotation_str[50] = {0};
                strncpy(rotation_str, row[rotation_col], sizeof(rotation_str)-1);

                /* Remove any quotes */
                pnp_remove_quotes(rotation_str);

                /* Try to parse rotation */
                if (1 != sscanf(rotation_str, "%lf", &pnpPartData.rotation)) {
                    pnpPartData.rotation = 0.0;
                }
            }

            /* Calculate rotation for package dimensions */
            gerb_transf_reset(tr_rot);
            gerb_transf_rotate(tr_rot, -DEG2RAD(pnpPartData.rotation));
            gerb_transf_apply(
                pnpPartData.pad_x - pnpPartData.mid_x, 
                pnpPartData.pad_y - pnpPartData.mid_y, 
                tr_rot, &tmp_x, &tmp_y
            );
        } else {
            printf("Debug: Skipping row with fewer than 3 columns\n");
            continue; /* Skip rows with too few columns */
        }

        /* Parse footprint shape */
        char* footprint_to_check = pnpPartData.footprint;
        int found_smd_pattern = 0;
        
        /* First check if we have a valid footprint from a properly identified footprint column */
        if (footprint_col >= 0 && 
            footprint_to_check != NULL && 
            strlen(footprint_to_check) > 0 && 
            strcmp(footprint_to_check, "(unknown)") != 0) {
            
            /* This comes from a column explicitly identified as 'footprint' or 'package' */
            printf("Debug: Using explicit footprint column: '%s'\n", footprint_to_check);
            
            /* Check if the footprint contains package dimensions */
            int direct_scan = sscanf(footprint_to_check, "%02d%02d", &i_length, &i_width);
            int prefix_scan = sscanf(footprint_to_check, "%*[^0-9]%02d%02d", &i_length, &i_width);
            
            printf("Debug: Footprint scan results - direct: %d, with prefix: %d\n", direct_scan, prefix_scan);
            
            if ((direct_scan == 2 || prefix_scan == 2) &&
                i_length >= 1 && i_length <= 25 && i_width >= 1 && i_width <= 12) {
                found_smd_pattern = 1;
                printf("Debug: Found valid SMD dimensions in footprint column\n");
            }
        }
        /* If no valid footprint pattern found, check comment/description column, but require "SMD" substring */
        else if (comment_col < ret && row[comment_col] != NULL && strstr(row[comment_col], "SMD")) {
            printf("Debug: Using comment/description column (index %d) for SMD footprint: '%s'\n", 
                   comment_col, row[comment_col]);
            footprint_to_check = row[comment_col];
            found_smd_pattern = 1;
        } else {
            printf("Debug: No valid footprint information found or field not identified as footprint\n");
        }
        
        printf("Debug: Analyzing final footprint: '%s'\n", footprint_to_check);
        
        /* Only proceed if we found a SMD pattern */
        if (found_smd_pattern) {
            /* If using description column, need to do dimension scanning again */
            if (footprint_to_check != pnpPartData.footprint) {
                int direct_scan = sscanf(footprint_to_check, "%02d%02d", &i_length, &i_width);
                int prefix_scan = sscanf(footprint_to_check, "%*[^0-9]%02d%02d", &i_length, &i_width);
                
                printf("Debug: Description scan results - direct: %d, with prefix: %d\n", direct_scan, prefix_scan);
                
                /* If scanning failed, don't treat as SMD */
                if (!((direct_scan == 2 || prefix_scan == 2) &&
                     i_length >= 1 && i_length <= 25 && i_width >= 1 && i_width <= 12)) {
                    found_smd_pattern = 0;
                }
            }
            
            if (found_smd_pattern) {
                
                /* Standard SMD packages (0603, 0805, etc) */
                printf("Debug: Detected SMD package: %02d%02d\n", i_length, i_width);
                pnpPartData.length = 0.01 * i_length;
                pnpPartData.width  = 0.01 * i_width;
                pnpPartData.shape  = PART_SHAPE_RECTANGLE;
            }
        }
        
        /* Handle mil-based footprints */
        else if (strstr(pnpPartData.footprint, "mil") && 
                sscanf(pnpPartData.footprint, "%d mil", &i_length) == 1) {

            /* Mil-based package dimensions */
            pnpPartData.length = i_length / 1000.0;
            pnpPartData.width = 0.1; /* Default 100 mil width */
            pnpPartData.shape = PART_SHAPE_RECTANGLE;

        } else {
            /* Default dimensions for unknown packages */
            pnpPartData.length = 0.1;  /* 100 mil square */
            pnpPartData.width = 0.1;   /* 100 mil square */
            pnpPartData.shape = PART_SHAPE_STD;
        }

        /* Add the part data to our array */
        g_array_append_val(pnpParseDataArray, pnpPartData);
        parsedLines += 1;
    }

    gerb_transf_free(tr_rot);

    /* Check if we parsed enough valid rows to consider this a PnP file */
    printf("Debug: Parsed %d lines out of %d total (%.1f%%)\n", 
           parsedLines, lineCounter, 100.0f * (float)parsedLines / (float)lineCounter);

    /* Consider the file valid if we found a valid header row OR have at least 3 parsed lines */
    if ((!foundValidDataRow && parsedLines < 3) || parsedLines == 0) {
        printf("Debug: Not enough valid data found (parsed %d of %d lines)\n", 
               parsedLines, lineCounter);
        g_array_free(pnpParseDataArray, TRUE);
        return NULL;
    }

    return pnpParseDataArray;
}

gboolean
pnp_check_and_parse_header(
    char* buf, int lineCounter, char* buf0, int* designator_col, int* footprint_col, 
    int* mid_x_col, int* mid_y_col, int* ref_x_col, int* ref_y_col, 
    int* pad_x_col, int* pad_y_col, int* layer_col, int* rotation_col, 
    int* comment_col, char* delimiter, gboolean* foundValidDataRow
) {
    /* Only check first 20 lines and stop if we've already found a header */
    if (lineCounter > 20 || *foundValidDataRow) {
        return FALSE;
    }

    /* Make a copy of the original buffer for parsing */
    char header_copy[MAXL + 2];
    strncpy(header_copy, buf, MAXL);
    header_copy[MAXL] = '\0';

    /* Skip comment characters at the beginning if present */
    char *header_start = header_copy;
    if (header_start[0] == '#' || header_start[0] == '*') {
        header_start++;
        /* Skip any whitespace after the comment character */
        while (*header_start && isspace(*header_start)) header_start++;
    }

    /* Check if this looks like a header line */
    if (*header_start) {

        printf("Debug: Check potential header row: %s\n", header_start);

        /* Clear temporary buffer */
        memset(buf0, 0, MAXL + 2);

        /* Auto-detect delimiter if not already set */
        if (*delimiter == '\0') {
            if (strstr(header_start, ",")) *delimiter = ',';
            else if (strstr(header_start, ";")) *delimiter = ';';
            else if (strstr(header_start, "|")) *delimiter = '|';
            else if (strstr(header_start, "\t")) *delimiter = '\t';
            else *delimiter = ','; /* Default to comma */

            printf("Debug: Auto-detected delimiter: '%c'\n", *delimiter);
        }

        /* Try to parse the header */
        if (pnp_parse_header_line(header_start, buf0, designator_col, footprint_col, 
                                mid_x_col, mid_y_col, ref_x_col, ref_y_col, 
                                pad_x_col, pad_y_col, layer_col, rotation_col, 
                                comment_col, delimiter)) {
            printf("Debug: Valid header row confirmed\n");
            *foundValidDataRow = TRUE;
            return TRUE;
        } else {
            printf("Debug: Not enough column names found, skipping as header\n");
        }
    }

    return FALSE;
}

/*	------------------------------------------------------------------
 *	pick_and_place_check_file_type
 *	------------------------------------------------------------------
 *	Description: Tries to parse the given file into a pick-and-place
 *		data set. If it fails to read any good rows, then returns
 *		FALSE, otherwise it returns TRUE.
 *	Notes:
 *	------------------------------------------------------------------
 */
gboolean
pick_and_place_check_file_type(gerb_file_t* fd, gboolean* returnFoundBinary) {
    char*    buf;
    int      len = 0;
    int      i;
    char*    letter;
    gboolean found_binary    = FALSE;
    gboolean found_G54       = FALSE;
    gboolean found_M0        = FALSE;
    gboolean found_M2        = FALSE;
    gboolean found_G2        = FALSE;
    gboolean found_ADD       = FALSE;
    gboolean found_comma     = FALSE;
    gboolean found_R         = FALSE;
    gboolean found_U         = FALSE;
    gboolean found_C         = FALSE;
    gboolean found_boardside = FALSE;

    buf = malloc(MAXL);
    if (buf == NULL)
        GERB_FATAL_ERROR("malloc buf failed in %s()", __FUNCTION__);

    while (fgets(buf, MAXL, fd->fd) != NULL) {
        len = strlen(buf);

        /* First look through the file for indications of its type */

        /* check for non-binary file */
        for (i = 0; i < len; i++) {
            if (!isprint((int)buf[i]) && (buf[i] != '\r') && (buf[i] != '\n') && (buf[i] != '\t')) {
                found_binary = TRUE;
            }
        }

        if (g_strstr_len(buf, len, "G54")) {
            found_G54 = TRUE;
        }
        if (g_strstr_len(buf, len, "M00")) {
            found_M0 = TRUE;
        }
        if (g_strstr_len(buf, len, "M02")) {
            found_M2 = TRUE;
        }
        if (g_strstr_len(buf, len, "G02")) {
            found_G2 = TRUE;
        }
        if (g_strstr_len(buf, len, "ADD")) {
            found_ADD = TRUE;
        }
        if (g_strstr_len(buf, len, ",")) {
            found_comma = TRUE;
        }
        /* Semicolon can be separator too */
        if (g_strstr_len(buf, len, ";")) {
            found_comma = TRUE;
        }

        /* Look for refdes -- This is dumb, but what else can we do? */
        if ((letter = g_strstr_len(buf, len, "R")) != NULL) {
            if (isdigit((int)letter[1])) { /* grab char after R */
                found_R = TRUE;
            }
        }
        if ((letter = g_strstr_len(buf, len, "C")) != NULL) {
            if (isdigit((int)letter[1])) { /* grab char after C */
                found_C = TRUE;
            }
        }
        if ((letter = g_strstr_len(buf, len, "U")) != NULL) {
            if (isdigit((int)letter[1])) { /* grab char after U */
                found_U = TRUE;
            }
        }

        /* Look for board side indicator since this is required
         * by many vendors */
        if (g_strstr_len(buf, len, "top")) {
            found_boardside = TRUE;
        }
        if (g_strstr_len(buf, len, "Top")) {
            found_boardside = TRUE;
        }
        if (g_strstr_len(buf, len, "TOP")) {
            found_boardside = TRUE;
        }
        /* Also look for evidence of "Layer" in header.... */
        if (g_strstr_len(buf, len, "ayer")) {
            found_boardside = TRUE;
        }
        if (g_strstr_len(buf, len, "AYER")) {
            found_boardside = TRUE;
        }
    }
    rewind(fd->fd);
    free(buf);

    /* Now form logical expression determining if this is a pick-place file */
    *returnFoundBinary = found_binary;
    if (found_G54)
        return FALSE;
    if (found_M0)
        return FALSE;
    if (found_M2)
        return FALSE;
    if (found_G2)
        return FALSE;
    if (found_ADD)
        return FALSE;
    if (found_comma && (found_R || found_C || found_U) && found_boardside)
        return TRUE;

    return FALSE;

} /* pick_and_place_check_file_type */

/*	------------------------------------------------------------------
 *	pick_and_place_convert_pnp_data_to_image
 *	------------------------------------------------------------------
 *	Description: Render a parsedPickAndPlaceData array into a gerb_image.
 *	Notes:
 *	------------------------------------------------------------------
 */
gerbv_image_t*
pick_and_place_convert_pnp_data_to_image(GArray* parsedPickAndPlaceData, gint boardSide) {
    gerbv_image_t*       image    = NULL;
    gerbv_net_t*         curr_net = NULL;
    gerbv_transf_t*      tr_rot   = gerb_transf_new();
    gerbv_drill_stats_t* stats; /* Eventually replace with pick_place_stats */
    gboolean             foundElement = FALSE;
    const double         draw_width   = 0.01;

    /* step through and make sure we have an element on the layer before
       we actually create a new image for it and fill it */
    for (guint i = 0; i < parsedPickAndPlaceData->len; i++) {
        PnpPartData partData = g_array_index(parsedPickAndPlaceData, PnpPartData, i);

        if ((boardSide == 0) && !((partData.layer[0] == 'b') || (partData.layer[0] == 'B')))
            continue;
        if ((boardSide == 1) && !((partData.layer[0] == 't') || (partData.layer[0] == 'T')))
            continue;

        foundElement = TRUE;
    }
    if (!foundElement)
        return NULL;

    image = gerbv_create_image(image, "Pick and Place (X-Y) File");
    if (image == NULL) {
        GERB_FATAL_ERROR("malloc image failed in %s()", __FUNCTION__);
    }

    image->format = g_new0(gerbv_format_t, 1);
    if (image->format == NULL) {
        GERB_FATAL_ERROR("malloc format failed in %s()", __FUNCTION__);
    }

    /* Separate top/bot layer type is needed for reload purpose */
    if (boardSide == 1)
        image->layertype = GERBV_LAYERTYPE_PICKANDPLACE_TOP;
    else
        image->layertype = GERBV_LAYERTYPE_PICKANDPLACE_BOT;

    stats = gerbv_drill_stats_new();
    if (stats == NULL)
        GERB_FATAL_ERROR("malloc pick_place_stats failed in %s()", __FUNCTION__);
    image->drill_stats = stats;

    curr_net        = image->netlist;
    curr_net->layer = image->layers;
    curr_net->state = image->states;
    pnp_reset_bbox(curr_net);
    image->info->min_x = HUGE_VAL;
    image->info->min_y = HUGE_VAL;
    image->info->max_x = -HUGE_VAL;
    image->info->max_y = -HUGE_VAL;

    image->aperture[0] = g_new0(gerbv_aperture_t, 1);
    assert(image->aperture[0] != NULL);
    image->aperture[0]->type           = GERBV_APTYPE_CIRCLE;
    image->aperture[0]->amacro         = NULL;
    image->aperture[0]->parameter[0]   = draw_width;
    image->aperture[0]->nuf_parameters = 1;

    for (guint i = 0; i < parsedPickAndPlaceData->len; i++) {
        PnpPartData partData = g_array_index(parsedPickAndPlaceData, PnpPartData, i);
        float       radius, labelOffset;

        curr_net        = pnp_new_net(curr_net);
        curr_net->layer = image->layers;
        curr_net->state = image->states;

        if ((partData.rotation > 89) && (partData.rotation < 91))
            labelOffset = fabs(partData.length / 2);
        else if ((partData.rotation > 179) && (partData.rotation < 181))
            labelOffset = fabs(partData.width / 2);
        else if ((partData.rotation > 269) && (partData.rotation < 271))
            labelOffset = fabs(partData.length / 2);
        else if ((partData.rotation > -91) && (partData.rotation < -89))
            labelOffset = fabs(partData.length / 2);
        else if ((partData.rotation > -181) && (partData.rotation < -179))
            labelOffset = fabs(partData.width / 2);
        else if ((partData.rotation > -271) && (partData.rotation < -269))
            labelOffset = fabs(partData.length / 2);
        else
            labelOffset = fabs(partData.width / 2);

        partData.rotation = DEG2RAD(partData.rotation);

        /* check if the entry is on the specified layer */
        if ((boardSide == 0) && !((partData.layer[0] == 'b') || (partData.layer[0] == 'B')))
            continue;
        if ((boardSide == 1) && !((partData.layer[0] == 't') || (partData.layer[0] == 'T')))
            continue;

        curr_net = pnp_new_net(curr_net);
        pnp_init_net(curr_net, image, partData.designator, GERBV_APERTURE_STATE_OFF, GERBV_INTERPOLATION_LINEARx1);

        /* First net of PNP is just a label holder, so calculate the lower left
         * location to line up above the element */
        curr_net->start_x = curr_net->stop_x = partData.mid_x;
        curr_net->start_y = curr_net->stop_y = partData.mid_y + labelOffset + draw_width;

        gerb_transf_reset(tr_rot);
        gerb_transf_shift(tr_rot, partData.mid_x, partData.mid_y);
        gerb_transf_rotate(tr_rot, -partData.rotation);

        if ((partData.shape == PART_SHAPE_RECTANGLE) || (partData.shape == PART_SHAPE_STD)) {
            // TODO: draw rectangle length x width taking into account rotation or pad x,y

            curr_net = pnp_new_net(curr_net);
            pnp_init_net(curr_net, image, partData.designator, GERBV_APERTURE_STATE_ON, GERBV_INTERPOLATION_LINEARx1);
            
            gerb_transf_apply(partData.length / 2, partData.width / 2, tr_rot, &curr_net->start_x, &curr_net->start_y);
            gerb_transf_apply(-partData.length / 2, partData.width / 2, tr_rot, &curr_net->stop_x, &curr_net->stop_y);

            /* TODO: write unifying function */

            curr_net = pnp_new_net(curr_net);
            pnp_init_net(curr_net, image, partData.designator, GERBV_APERTURE_STATE_ON, GERBV_INTERPOLATION_LINEARx1);

            gerb_transf_apply(-partData.length / 2, partData.width / 2, tr_rot, &curr_net->start_x, &curr_net->start_y);
            gerb_transf_apply(-partData.length / 2, -partData.width / 2, tr_rot, &curr_net->stop_x, &curr_net->stop_y);

            curr_net = pnp_new_net(curr_net);
            pnp_init_net(curr_net, image, partData.designator, GERBV_APERTURE_STATE_ON, GERBV_INTERPOLATION_LINEARx1);

            gerb_transf_apply(
                -partData.length / 2, -partData.width / 2, tr_rot, &curr_net->start_x, &curr_net->start_y
            );
            gerb_transf_apply(partData.length / 2, -partData.width / 2, tr_rot, &curr_net->stop_x, &curr_net->stop_y);

            curr_net = pnp_new_net(curr_net);
            pnp_init_net(curr_net, image, partData.designator, GERBV_APERTURE_STATE_ON, GERBV_INTERPOLATION_LINEARx1);

            gerb_transf_apply(partData.length / 2, -partData.width / 2, tr_rot, &curr_net->start_x, &curr_net->start_y);
            gerb_transf_apply(partData.length / 2, partData.width / 2, tr_rot, &curr_net->stop_x, &curr_net->stop_y);

            curr_net = pnp_new_net(curr_net);
            pnp_init_net(curr_net, image, partData.designator, GERBV_APERTURE_STATE_ON, GERBV_INTERPOLATION_LINEARx1);

            if (partData.shape == PART_SHAPE_RECTANGLE) {
                gerb_transf_apply(
                    partData.length / 4, -partData.width / 2, tr_rot, &curr_net->start_x, &curr_net->start_y
                );
                gerb_transf_apply(
                    partData.length / 4, partData.width / 2, tr_rot, &curr_net->stop_x, &curr_net->stop_y
                );
            } else {
                gerb_transf_apply(
                    partData.length / 4, partData.width / 2, tr_rot, &curr_net->start_x, &curr_net->start_y
                );
                gerb_transf_apply(
                    partData.length / 4, partData.width / 4, tr_rot, &curr_net->stop_x, &curr_net->stop_y
                );

                curr_net = pnp_new_net(curr_net);
                pnp_init_net(
                    curr_net, image, partData.designator, GERBV_APERTURE_STATE_ON, GERBV_INTERPOLATION_LINEARx1
                );

                gerb_transf_apply(
                    partData.length / 2, partData.width / 4, tr_rot, &curr_net->start_x, &curr_net->start_y
                );
                gerb_transf_apply(
                    partData.length / 4, partData.width / 4, tr_rot, &curr_net->stop_x, &curr_net->stop_y
                );
            }

            /* calculate a rough radius for the min/max screen calcs later */
            radius = MAX(partData.length / 2, partData.width / 2);
        } else {
            gdouble tmp_x, tmp_y;

            pnp_init_net(curr_net, image, partData.designator, GERBV_APERTURE_STATE_ON, GERBV_INTERPOLATION_LINEARx1);

            curr_net->start_x = partData.mid_x;
            curr_net->start_y = partData.mid_y;
            gerb_transf_apply(partData.pad_x - partData.mid_x, partData.pad_y - partData.mid_y, tr_rot, &tmp_x, &tmp_y);

            curr_net->stop_x = tmp_x;
            curr_net->stop_y = tmp_y;

            curr_net = pnp_new_net(curr_net);
            pnp_init_net(
                curr_net, image, partData.designator, GERBV_APERTURE_STATE_ON, GERBV_INTERPOLATION_CW_CIRCULAR
            );

            curr_net->start_x = partData.mid_x;
            curr_net->start_y = partData.mid_y;
            curr_net->stop_x  = partData.pad_x;
            curr_net->stop_y  = partData.pad_y;

            curr_net->cirseg         = g_new0(gerbv_cirseg_t, 1);
            curr_net->cirseg->angle1 = 0.0;
            curr_net->cirseg->angle2 = 360.0;
            curr_net->cirseg->cp_x   = partData.mid_x;
            curr_net->cirseg->cp_y   = partData.mid_y;
            radius                   = hypot(partData.pad_x - partData.mid_x, partData.pad_y - partData.mid_y);
            if (radius < 0.001)
                radius = 0.1;
            curr_net->cirseg->width  = 2 * radius; /* fabs(pad_x-mid_x) */
            curr_net->cirseg->height = 2 * radius;
        }

        /*
         * update min and max numbers so the screen zoom-to-fit
         *function will work
         */
        image->info->min_x = MIN(image->info->min_x, (partData.mid_x - radius - 0.02));
        image->info->min_y = MIN(image->info->min_y, (partData.mid_y - radius - 0.02));
        image->info->max_x = MAX(image->info->max_x, (partData.mid_x + radius + 0.02));
        image->info->max_y = MAX(image->info->max_y, (partData.mid_y + radius + 0.02));
    }
    curr_net->next = NULL;

    gerb_transf_free(tr_rot);
    return image;
} /* pick_and_place_convert_pnp_data_to_image */

/*	------------------------------------------------------------------
 *	pick_and_place_parse_file_to_images
 *	------------------------------------------------------------------
 *	Description: Renders a pick and place file to a gerb_image.
 *	If image pointer is not NULL, then corresponding image will not be
 *	populated.
 *	Notes: The file format should already be verified before calling
 *       this function, since it does very little sanity checking itself.
 *	------------------------------------------------------------------
 */
void
pick_and_place_parse_file_to_images(gerb_file_t* fd, gerbv_image_t** topImage, gerbv_image_t** bottomImage) {
    GArray* parsedPickAndPlaceData = pick_and_place_parse_file(fd);

    if (parsedPickAndPlaceData != NULL) {
        /* Non NULL pointer is used as "not to reload" mark */
        if (*bottomImage == NULL)
            *bottomImage = pick_and_place_convert_pnp_data_to_image(parsedPickAndPlaceData, 0);

        if (*topImage == NULL)
            *topImage = pick_and_place_convert_pnp_data_to_image(parsedPickAndPlaceData, 1);

        g_array_free(parsedPickAndPlaceData, TRUE);
    }
} /* pick_and_place_parse_file_to_images */

static gerbv_net_t*
pnp_new_net(gerbv_net_t* net) {
    gerbv_net_t* n;
    net->next = g_new0(gerbv_net_t, 1);
    n         = net->next;
    assert(n != NULL);

    pnp_reset_bbox(n);

    return n;
}

static void
pnp_reset_bbox(gerbv_net_t* net) {
    net->boundingBox.left   = -HUGE_VAL;
    net->boundingBox.right  = HUGE_VAL;
    net->boundingBox.bottom = -HUGE_VAL;
    net->boundingBox.top    = HUGE_VAL;
}

static void
pnp_init_net(
    gerbv_net_t* net, gerbv_image_t* image, const char* label, gerbv_aperture_state_t apert_state,
    gerbv_interpolation_t interpol
) {
    net->aperture       = 0;
    net->aperture_state = apert_state;
    net->interpolation  = interpol;
    net->layer          = image->layers;
    net->state          = image->states;

    if (strlen(label) > 0) {
        net->label = g_string_new(label);
    }
}

/**
 * Parses a header line to identify column positions.
 * Handles both standard headers and commented headers.
 * 
 * @param buf Input line buffer
 * @param buf0 Temporary buffer
 * @param designator_col Pointer to designator column index
 * @param footprint_col Pointer to footprint column index
 * @param mid_x_col Pointer to X coordinate column index
 * @param mid_y_col Pointer to Y coordinate column index
 * @param ref_x_col Pointer to reference X column index
 * @param ref_y_col Pointer to reference Y column index
 * @param pad_x_col Pointer to pad X column index
 * @param pad_y_col Pointer to pad Y column index
 * @param layer_col Pointer to layer column index
 * @param rotation_col Pointer to rotation column index
 * @param comment_col Pointer to comment column index
 * @param delimiter Pointer to delimiter character
 * 
 * @return TRUE if header was successfully parsed, FALSE otherwise
 */
/** 
 * Custom function to parse a header line that is comma-delimited.
 * We need this because the csv_row_parse function has issues with our header format.
 * This function handles leading/trailing whitespace and returns an array of field pointers.
 */
static int 
custom_parse_comma_header(char* input, char* result_buffer, char* fields[], int max_fields) {
    char* src = input;
    char* dest = result_buffer;
    int field_count = 0;
    gboolean in_quote = FALSE;
    
    /* Skip leading whitespace */
    while (*src && isspace(*src)) src++;
    
    /* Set the first field */
    fields[field_count++] = dest;
    
    /* Process each character */
    while (*src && field_count < max_fields) {
        if (*src == '"') {
            /* Handle quotes - toggle quote state but don't include them in output */
            in_quote = !in_quote;
            src++;  /* Skip the quote character */
        } else if (*src == ',' && !in_quote) {
            /* End of field (only if not inside quotes) */
            *dest++ = '\0';  /* Terminate the current field */
            
            /* Skip the delimiter and any whitespace after it */
            src++;
            while (*src && isspace(*src) && *src != '"') src++;
            
            /* Start a new field if we haven't reached the max */
            if (field_count < max_fields) {
                fields[field_count++] = dest;
            }
        } else {
            /* Copy the character */
            *dest++ = *src++;
        }
    }
    
    /* Terminate the last field */
    *dest = '\0';
    
    /* Trim trailing whitespace from all fields */
    for (int i = 0; i < field_count; i++) {
        int len = strlen(fields[i]);
        while (len > 0 && isspace(fields[i][len-1])) {
            fields[i][--len] = '\0';
        }
    }
    
    printf("Debug: Custom parser found %d fields\n", field_count);
    
    return field_count;
}

static gboolean
pnp_parse_header_line(
    char* buf, char* buf0, int* designator_col, int* footprint_col, int* mid_x_col, int* mid_y_col,
    int* ref_x_col, int* ref_y_col, int* pad_x_col, int* pad_y_col, int* layer_col, int* rotation_col,
    int* comment_col, char* delimiter
) {
    char header_buf[MAXL + 2];
    char* start = buf;
    int header_ret = 0;
    char* header_row[12] = {NULL}; /* Initialize all pointers to NULL */
    
    /* Initialize header buffer with zeros to ensure clean parsing */
    memset(header_buf, 0, MAXL + 2);
    memset(buf0, 0, MAXL + 2);
    
    /* Handle commented headers - skip # and any whitespace */
    if (buf[0] == '#') {
        start = buf + 1;
        while (*start && isspace(*start)) start++;
    }
    
    /* For safety, limit copying to visible ASCII characters only */
    int i = 0, j = 0;
    while (start[i] && j < MAXL - 1) {
        /* Skip control characters and extended ASCII */
        if (start[i] >= 32 && start[i] <= 126) {
            header_buf[j++] = start[i];
        }
        i++;
    }
    header_buf[j] = '\0'; /* Ensure null termination */
    
    /* Check if the buffer is empty after sanitization */
    if (header_buf[0] == '\0') {
        printf("Debug: Header line is empty after sanitization\n");
        return FALSE;
    }
    
    /* Use comma as default delimiter but check for others */
    if (*delimiter == '\0') {
        /* Count delimiters */
        int comma_count = 0, tab_count = 0, semicolon_count = 0, pipe_count = 0;
        
        for (i = 0; header_buf[i]; i++) {
            if (header_buf[i] == ',') comma_count++;
            if (header_buf[i] == '\t') tab_count++;
            if (header_buf[i] == ';') semicolon_count++;
            if (header_buf[i] == '|') pipe_count++;
        }
        
        printf("Debug: Found delimiters in header - commas: %d, tabs: %d, semicolons: %d, pipes: %d\n",
               comma_count, tab_count, semicolon_count, pipe_count);
        
        /* Choose the most frequent delimiter */
        *delimiter = ',';  /* Default to comma */
        int max_count = comma_count;
        
        if (tab_count > max_count) {
            max_count = tab_count;
            *delimiter = '\t';
        }
        if (semicolon_count > max_count) {
            max_count = semicolon_count;
            *delimiter = ';';
        }
        if (pipe_count > max_count) {
            *delimiter = '|';
        }
        
        printf("Debug: Selected delimiter: '%c'\n", *delimiter);
    }
    
    /* For comma delimiter, use our custom parser */
    if (*delimiter == ',') {
        printf("Debug: Using custom parser for comma-delimited header: '%s'\n", header_buf);
        header_ret = custom_parse_comma_header(header_buf, buf0, header_row, 11);
    } else {
        /* For other delimiters, try the standard parser but with safety checks */
        printf("Debug: Using standard parser with delimiter '%c' for: '%s'\n", *delimiter, header_buf);
        header_ret = csv_row_parse(header_buf, strlen(header_buf), buf0, MAXL, header_row, 11, (int)*delimiter, CSV_QUOTES);
        
        /* Handle error from csv_row_parse */
        if (header_ret <= 0) {
            printf("Debug: Standard parser failed, trying custom comma parser as fallback\n");
            header_ret = custom_parse_comma_header(header_buf, buf0, header_row, 11);
        }
    }
    
    /* Print all returned fields for debug */
    printf("Debug: Header parse returned %d fields\n", header_ret);
    for (i = 0; i < header_ret && i < 11; i++) {
        if (header_row[i] != NULL) {
            printf("Debug: Field[%d] = '%s'\n", i, header_row[i]);
        } else {
            printf("Debug: Field[%d] = NULL\n", i);
        }
    }
    
    /* Verify this is actually a header row by counting known column names */
    int known_column_count = 0;
    
    if (header_ret <= 0) {
        return FALSE;
    }
    
    /* Try to identify column positions */
    for (i = 0; i < header_ret && i < 11; i++) {
        if (header_row[i] == NULL || header_row[i][0] == '\0') {
            printf("Debug: Header[%d] is empty or NULL\n", i);
            continue;
        }
        
        printf("Debug: Processing header[%d] = '%s'\n", i, header_row[i]);
        
        /* Convert to lowercase for case-insensitive matching */
        char temp_buf[50];
        if (strlen(header_row[i]) >= sizeof(temp_buf)) {
            /* Too long, only copy what we can fit */
            strncpy(temp_buf, header_row[i], sizeof(temp_buf)-1);
            temp_buf[sizeof(temp_buf)-1] = '\0';
        } else {
            strcpy(temp_buf, header_row[i]);
        }
        
        /* Trim leading whitespace */
        char* p = temp_buf;
        while (*p && isspace(*p)) p++;
        
        /* If necessary, shift the string left */
        if (p > temp_buf) {
            memmove(temp_buf, p, strlen(p) + 1);
        }
        
        /* Convert to lowercase */
        for (p = temp_buf; *p; ++p) *p = tolower(*p);
        
        /* Check for known column names - being more strict now */
        if ((strcmp(temp_buf, "refdes") == 0) || 
            (strcmp(temp_buf, "ref") == 0) || 
            (strncmp(temp_buf, "des", 3) == 0 && !strstr(temp_buf, "description"))) {
            *designator_col = i;
            printf("Debug: Found designator column at index %d\n", i);
            known_column_count++;
        } else if (strcmp(temp_buf, "description") == 0) {
            /* Map 'Description' to the comment field */
            *comment_col = i;
            printf("Debug: Found description column at index %d (mapping to comment)\n", i);
            known_column_count++;
        } else if (strstr(temp_buf, "foot") || strstr(temp_buf, "pack")) {
            /* Only use foot(print) or pack(age) fields as footprint */
            *footprint_col = i;
            printf("Debug: Found footprint/package column at index %d\n", i);
            known_column_count++;
        } else if (strcmp(temp_buf, "value") == 0) {
            /* Sometimes Value contains the part value, not footprint info */
            /* Do NOT treat value as footprint by default */
            printf("Debug: Found value column at index %d (not using as footprint by default)\n", i);
        } else if ((strstr(temp_buf, "mid") && strstr(temp_buf, "x")) || 
                   strcmp(temp_buf, "x") == 0) {
            *mid_x_col = i;
            printf("Debug: Found mid_x column at index %d\n", i);
            known_column_count++;
        } else if ((strstr(temp_buf, "mid") && strstr(temp_buf, "y")) || 
                   strcmp(temp_buf, "y") == 0) {
            *mid_y_col = i;
            printf("Debug: Found mid_y column at index %d\n", i);
            known_column_count++;
        } else if (strstr(temp_buf, "rot")) {
            *rotation_col = i;
            printf("Debug: Found rotation column at index %d\n", i);
            known_column_count++;
        } else if (strstr(temp_buf, "layer") || strstr(temp_buf, "side") || 
                  strstr(temp_buf, "top") || strstr(temp_buf, "bottom") ||
                  strstr(temp_buf, "top/bottom")) {
            *layer_col = i;
            printf("Debug: Found layer column at index %d\n", i);
            known_column_count++;
        } else if (strstr(temp_buf, "comm") || strstr(temp_buf, "val")) {
            *comment_col = i;
            printf("Debug: Found comment column at index %d\n", i);
            known_column_count++;
        }
    }
    
    printf("Debug: Found %d known column names\n", known_column_count);
    
    /* Special case for PCB XY format with specific header format */
    if (strstr(buf, "RefDes") && strstr(buf, "X") && strstr(buf, "Y") && 
        strstr(buf, "rotation") && known_column_count < 3) {
        
        printf("Debug: Recognized PCB XY format header despite low column count\n");
        
        /* Set columns based on pattern matching in the header */
        for (i = 0; i < header_ret; i++) {
            if (header_row[i] == NULL) continue;
            
            if (strstr(header_row[i], "RefDes"))
                *designator_col = i;
            else if (strstr(header_row[i], "Footprint") || strstr(header_row[i], "Package"))
                *footprint_col = i;
            else if (strcmp(header_row[i], "X") == 0 || strcmp(header_row[i], " X") == 0)
                *mid_x_col = i;
            else if (strcmp(header_row[i], "Y") == 0 || strcmp(header_row[i], " Y") == 0)
                *mid_y_col = i;
            else if (strstr(header_row[i], "rotation"))
                *rotation_col = i;
            else if (strstr(header_row[i], "top") || strstr(header_row[i], "bottom"))
                *layer_col = i;
        }
        
        printf("Debug: Using PCB XY pattern-matched column map: RefDes=%d Value=%d X=%d Y=%d Rot=%d Layer=%d\n",
               *designator_col, *footprint_col, *mid_x_col, *mid_y_col, *rotation_col, *layer_col);
               
        return TRUE;
    }
    
    /* Use the specific PCB XY format pattern */
    if (strstr(buf, "# RefDes") && 
        (strstr(buf, ", Description") || strstr(buf, " Description")) && 
        (strstr(buf, ", Value") || strstr(buf, " Value")) && 
        ((strstr(buf, ", X") || strstr(buf, " X")) && 
         (strstr(buf, ", Y") || strstr(buf, " Y"))) && 
        (strstr(buf, ", rotation") || strstr(buf, " rotation")) && 
        (strstr(buf, ", top/bottom") || strstr(buf, " top/bottom"))) {
        
        printf("Debug: Matched exact PCB XY header format\n");
        
        /* Set columns to standard PCB XY format */
        *designator_col = 0;
        *footprint_col = 2;  /* Value column typically used for footprint */
        *mid_x_col = 3;
        *mid_y_col = 4;
        *rotation_col = 5;
        *layer_col = 6;
        
        printf("Debug: Using standard PCB XY column map: RefDes=%d Value=%d X=%d Y=%d Rot=%d Layer=%d\n",
               *designator_col, *footprint_col, *mid_x_col, *mid_y_col, *rotation_col, *layer_col);
        
        return TRUE;
    }
    
    /* Only consider it a valid header if we found at least 3 known column names */
    return (known_column_count >= 3);
}
