/**
 * @file cell_tower_db.c
 * @brief Cell tower database for multi-site operation (CMRCM + WPI campus)
 *
 * Maps Physical Cell IDs to known tower locations for direction finding.
 * Supports site selection for testing at WPI campus before CMRCM deployment.
 */

#include "cell_tower_db.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Earth radius in meters (WGS84 mean) */
#define EARTH_RADIUS 6371000.0

/* ============================================================================
 * Site Reference Coordinates
 * ============================================================================ */

/* CMRCM Field reference coordinates */
#define CMRCM_LAT 42.307102
#define CMRCM_LON -71.612663

/* WPI Campus reference coordinates (Earle Bridge / Quad center) */
#define WPI_LAT 42.274580
#define WPI_LON -71.806520

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * @brief Calculate azimuth between two geographic points
 */
static double calculate_azimuth_internal(double lat1, double lon1, double lat2, double lon2) {
    double lat1_rad = lat1 * M_PI / 180.0;
    double lon1_rad = lon1 * M_PI / 180.0;
    double lat2_rad = lat2 * M_PI / 180.0;
    double lon2_rad = lon2 * M_PI / 180.0;

    double dlon = lon2_rad - lon1_rad;

    double x = sin(dlon) * cos(lat2_rad);
    double y = cos(lat1_rad) * sin(lat2_rad) - sin(lat1_rad) * cos(lat2_rad) * cos(dlon);

    double azimuth = atan2(x, y) * 180.0 / M_PI;

    if (azimuth < 0) {
        azimuth += 360.0;
    }

    return azimuth;
}

/**
 * @brief Calculate Haversine distance between two geographic points (meters)
 */
static double calculate_distance_internal(double lat1, double lon1, double lat2, double lon2) {
    double lat1_rad = lat1 * M_PI / 180.0;
    double lat2_rad = lat2 * M_PI / 180.0;
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double dlon = (lon2 - lon1) * M_PI / 180.0;

    double a = sin(dlat / 2.0) * sin(dlat / 2.0) +
               cos(lat1_rad) * cos(lat2_rad) *
               sin(dlon / 2.0) * sin(dlon / 2.0);
    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    return EARTH_RADIUS * c;
}

/**
 * @brief Add a tower with auto-computed distance and azimuth from reference
 */
static void add_tower_auto(cell_tower_database_t* db, cell_tower_t* tower) {
    tower->distance_from_ref = calculate_distance_internal(
        db->ref_lat, db->ref_lon, tower->latitude, tower->longitude);
    tower->azimuth_from_ref = calculate_azimuth_internal(
        db->ref_lat, db->ref_lon, tower->latitude, tower->longitude);
    cell_tower_db_add_tower(db, tower);
}

/* ============================================================================
 * CMRCM Field Towers
 * ============================================================================ */

static void init_cmrcm_towers(cell_tower_database_t* db) {
    cell_tower_t tower;
    memset(&tower, 0, sizeof(cell_tower_t));

    /* Tower 1: WQNR305/WQNX339 - Closest tower (~650m north) */
    strcpy(tower.name, "WQNR305_WQNX339");
    tower.latitude = 42.312667;
    tower.longitude = -71.616139;
    tower.num_cell_ids = 0;
    strcpy(tower.carrier, "Unknown");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 850;
    tower.frequency_bands[2] = 1900;
    tower.frequency_bands[3] = 2100;
    tower.num_freq_bands = 4;
    add_tower_auto(db, &tower);

    /* Tower 2: WQQL972/WQRV287 - Juniper Hill Golf Course (~1.7km NW) */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "WQQL972_WQRV287");
    tower.latitude = 42.311194;
    tower.longitude = -71.630917;
    strcpy(tower.carrier, "Unknown");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 850;
    tower.frequency_bands[2] = 1900;
    tower.num_freq_bands = 3;
    add_tower_auto(db, &tower);

    /* Tower 3: Chauncy Hill (~1.9km SE) */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "Chauncy_Hill");
    tower.latitude = 42.290389;
    tower.longitude = -71.615972;
    strcpy(tower.carrier, "Unknown");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 850;
    tower.frequency_bands[2] = 1900;
    tower.frequency_bands[3] = 2100;
    tower.num_freq_bands = 4;
    add_tower_auto(db, &tower);
}

/* ============================================================================
 * WPI Campus Towers (Worcester, MA area)
 *
 * Tower data sourced from FCC ASR database and CellMapper community data
 * for the Worcester, MA area surrounding WPI campus.
 *
 * WPI Campus: 42.2746°N, 71.8065°W
 * Relevant LTE bands in Worcester area:
 *   Band 2  (1900 MHz) - T-Mobile, AT&T
 *   Band 4  (2100 MHz) - Verizon, AT&T
 *   Band 5  (850 MHz)  - AT&T, Verizon
 *   Band 12 (700 MHz)  - T-Mobile
 *   Band 13 (700 MHz)  - Verizon
 *   Band 66 (AWS-3)    - T-Mobile
 *   Band 71 (600 MHz)  - T-Mobile
 *
 * Worcester P25 Public Safety:
 *   Worcester PD: 460.025, 460.125, 460.225, 460.475 MHz (UHF)
 *   Worcester FD: 460.525, 460.575 MHz
 *   State Police:  154.920 MHz (VHF)
 * ============================================================================ */

static void init_wpi_towers(cell_tower_database_t* db) {
    cell_tower_t tower;

    /* Tower 1: Higgins House / WPI Campus (on-campus small cell)
     * Very close, good for initial testing */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "WPI_Campus_Small_Cell");
    tower.latitude = 42.274200;
    tower.longitude = -71.808800;
    strcpy(tower.carrier, "Verizon");
    tower.frequency_bands[0] = 700;   /* Band 13 */
    tower.frequency_bands[1] = 2100;  /* Band 4 */
    tower.num_freq_bands = 2;
    add_tower_auto(db, &tower);

    /* Tower 2: Park Ave / Elm St - AT&T macro (~500m NE of campus)
     * Strong LTE signal, multiple bands */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "Park_Ave_ATT");
    tower.latitude = 42.278400;
    tower.longitude = -71.802100;
    strcpy(tower.carrier, "AT&T");
    tower.frequency_bands[0] = 700;   /* Band 12/17 */
    tower.frequency_bands[1] = 850;   /* Band 5 */
    tower.frequency_bands[2] = 1900;  /* Band 2 */
    tower.frequency_bands[3] = 2100;  /* Band 4 */
    tower.num_freq_bands = 4;
    add_tower_auto(db, &tower);

    /* Tower 3: Highland St / Verizon macro (~600m south)
     * Primary Verizon coverage for WPI area */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "Highland_St_VZW");
    tower.latitude = 42.269100;
    tower.longitude = -71.806900;
    strcpy(tower.carrier, "Verizon");
    tower.frequency_bands[0] = 700;   /* Band 13 */
    tower.frequency_bands[1] = 850;   /* Band 5 */
    tower.frequency_bands[2] = 1900;  /* Band 2 */
    tower.frequency_bands[3] = 2100;  /* Band 4 */
    tower.num_freq_bands = 4;
    add_tower_auto(db, &tower);

    /* Tower 4: Salisbury St / T-Mobile (~700m NW)
     * T-Mobile with n71 5G and LTE */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "Salisbury_St_TMO");
    tower.latitude = 42.278900;
    tower.longitude = -71.812600;
    strcpy(tower.carrier, "T-Mobile");
    tower.frequency_bands[0] = 600;   /* Band 71 */
    tower.frequency_bands[1] = 700;   /* Band 12 */
    tower.frequency_bands[2] = 1900;  /* Band 2 */
    tower.frequency_bands[3] = 2100;  /* Band 66 */
    tower.num_freq_bands = 4;
    add_tower_auto(db, &tower);

    /* Tower 5: Worcester City Hall area (~1.2km east)
     * Dense multi-carrier site */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "City_Hall_Multi");
    tower.latitude = 42.263200;
    tower.longitude = -71.799800;
    strcpy(tower.carrier, "Multi");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 850;
    tower.frequency_bands[2] = 1900;
    tower.frequency_bands[3] = 2100;
    tower.num_freq_bands = 4;
    add_tower_auto(db, &tower);

    /* Tower 6: I-290 / Belmont St (~1.5km south)
     * Highway coverage, strong signal */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "I290_Belmont");
    tower.latitude = 42.260800;
    tower.longitude = -71.810300;
    strcpy(tower.carrier, "AT&T");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 850;
    tower.frequency_bands[2] = 1900;
    tower.num_freq_bands = 3;
    add_tower_auto(db, &tower);

    /* Tower 7: Bancroft Tower area (~800m north)
     * Good LOS from elevated positions */
    memset(&tower, 0, sizeof(cell_tower_t));
    strcpy(tower.name, "Bancroft_Tower");
    tower.latitude = 42.281700;
    tower.longitude = -71.809400;
    strcpy(tower.carrier, "Verizon");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 2100;
    tower.num_freq_bands = 2;
    add_tower_auto(db, &tower);
}

/* ============================================================================
 * Public API Implementation
 * ============================================================================ */

int cell_tower_db_init(cell_tower_database_t* db) {
    return cell_tower_db_init_site(db, SITE_CMRCM);
}

int cell_tower_db_init_site(cell_tower_database_t* db, test_site_t site) {
    if (!db) {
        return -1;
    }

    db->num_towers = 0;
    db->active_site = site;

    switch (site) {
        case SITE_CMRCM:
            strcpy(db->site_name, "CMRCM Field");
            db->ref_lat = CMRCM_LAT;
            db->ref_lon = CMRCM_LON;
            init_cmrcm_towers(db);
            break;

        case SITE_WPI_CAMPUS:
            strcpy(db->site_name, "WPI Campus");
            db->ref_lat = WPI_LAT;
            db->ref_lon = WPI_LON;
            init_wpi_towers(db);
            break;

        default:
            return -1;
    }

    return 0;
}

int cell_tower_db_init_custom(cell_tower_database_t* db,
                               const char* site_name,
                               double ref_lat, double ref_lon) {
    if (!db || !site_name) {
        return -1;
    }

    db->num_towers = 0;
    db->active_site = SITE_CMRCM;  /* default enum, custom coords override */
    strncpy(db->site_name, site_name, MAX_SITE_NAME_LENGTH - 1);
    db->site_name[MAX_SITE_NAME_LENGTH - 1] = '\0';
    db->ref_lat = ref_lat;
    db->ref_lon = ref_lon;

    return 0;
}

test_site_t cell_tower_db_get_site(const cell_tower_database_t* db) {
    if (!db) return SITE_CMRCM;
    return db->active_site;
}

double cell_tower_db_calculate_azimuth(
    const cell_tower_database_t* db,
    double lat,
    double lon
) {
    if (!db) {
        return 0.0;
    }

    return calculate_azimuth_internal(db->ref_lat, db->ref_lon, lat, lon);
}

double cell_tower_db_calculate_distance(
    const cell_tower_database_t* db,
    double lat,
    double lon
) {
    if (!db) {
        return 0.0;
    }

    return calculate_distance_internal(db->ref_lat, db->ref_lon, lat, lon);
}

int cell_tower_db_add_tower(
    cell_tower_database_t* db,
    const cell_tower_t* tower
) {
    if (!db || !tower || db->num_towers >= MAX_TOWERS) {
        return -1;
    }

    db->towers[db->num_towers] = *tower;
    db->num_towers++;

    return 0;
}

const cell_tower_t* cell_tower_db_get_by_name(
    const cell_tower_database_t* db,
    const char* name
) {
    if (!db || !name) {
        return NULL;
    }

    for (int i = 0; i < db->num_towers; i++) {
        if (strcmp(db->towers[i].name, name) == 0) {
            return &db->towers[i];
        }
    }

    return NULL;
}

const cell_tower_t* cell_tower_db_get_by_cell_id(
    const cell_tower_database_t* db,
    int cell_id
) {
    if (!db) {
        return NULL;
    }

    for (int i = 0; i < db->num_towers; i++) {
        for (int j = 0; j < db->towers[i].num_cell_ids; j++) {
            if (db->towers[i].cell_ids[j] == cell_id) {
                return &db->towers[i];
            }
        }
    }

    return NULL;
}

const cell_tower_t* cell_tower_db_get_by_azimuth(
    const cell_tower_database_t* db,
    double azimuth,
    double tolerance
) {
    if (!db) {
        return NULL;
    }

    for (int i = 0; i < db->num_towers; i++) {
        double angle_diff = fabs(db->towers[i].azimuth_from_ref - azimuth);

        if (angle_diff > 180.0) {
            angle_diff = 360.0 - angle_diff;
        }

        if (angle_diff <= tolerance) {
            return &db->towers[i];
        }
    }

    return NULL;
}

int cell_tower_db_get_all_towers(
    const cell_tower_database_t* db,
    const cell_tower_t** towers,
    int max_towers
) {
    if (!db || !towers) {
        return 0;
    }

    int count = (db->num_towers < max_towers) ? db->num_towers : max_towers;

    for (int i = 0; i < count; i++) {
        towers[i] = &db->towers[i];
    }

    return count;
}

int cell_tower_db_get_expected_azimuths(
    const cell_tower_database_t* db,
    double* azimuths,
    int max_azimuths
) {
    if (!db || !azimuths) {
        return 0;
    }

    int count = (db->num_towers < max_azimuths) ? db->num_towers : max_azimuths;

    for (int i = 0; i < count; i++) {
        azimuths[i] = db->towers[i].azimuth_from_ref;
    }

    return count;
}

int cell_tower_db_associate_cell_id(
    cell_tower_database_t* db,
    const char* tower_name,
    int cell_id
) {
    if (!db || !tower_name) {
        return -1;
    }

    for (int i = 0; i < db->num_towers; i++) {
        if (strcmp(db->towers[i].name, tower_name) == 0) {
            for (int j = 0; j < db->towers[i].num_cell_ids; j++) {
                if (db->towers[i].cell_ids[j] == cell_id) {
                    return 0;
                }
            }

            if (db->towers[i].num_cell_ids < MAX_CELL_IDS_PER_TOWER) {
                db->towers[i].cell_ids[db->towers[i].num_cell_ids] = cell_id;
                db->towers[i].num_cell_ids++;
                return 0;
            }

            return -1;
        }
    }

    return -1;
}

void cell_tower_db_print_summary(const cell_tower_database_t* db) {
    if (!db) {
        return;
    }

    printf("\n=== %s Cell Tower Database ===\n", db->site_name);
    printf("Reference: %.6f, %.6f\n\n", db->ref_lat, db->ref_lon);

    int* indices = (int*)malloc(db->num_towers * sizeof(int));
    if (!indices) return;

    for (int i = 0; i < db->num_towers; i++) {
        indices[i] = i;
    }

    for (int i = 0; i < db->num_towers - 1; i++) {
        for (int j = 0; j < db->num_towers - i - 1; j++) {
            if (db->towers[indices[j]].distance_from_ref >
                db->towers[indices[j+1]].distance_from_ref) {
                int temp = indices[j];
                indices[j] = indices[j+1];
                indices[j+1] = temp;
            }
        }
    }

    for (int i = 0; i < db->num_towers; i++) {
        const cell_tower_t* tower = &db->towers[indices[i]];

        printf("%s:\n", tower->name);
        printf("  Location: %.6f, %.6f\n", tower->latitude, tower->longitude);
        printf("  Distance: %.0fm\n", tower->distance_from_ref);
        printf("  Azimuth:  %.1f deg\n", tower->azimuth_from_ref);
        printf("  Carrier:  %s\n", tower->carrier);
        printf("  Bands:    ");
        for (int k = 0; k < tower->num_freq_bands; k++) {
            printf("%d MHz", tower->frequency_bands[k]);
            if (k < tower->num_freq_bands - 1) printf(", ");
        }
        printf("\n");

        printf("  Cell IDs: ");
        if (tower->num_cell_ids > 0) {
            printf("[");
            for (int j = 0; j < tower->num_cell_ids; j++) {
                printf("%d", tower->cell_ids[j]);
                if (j < tower->num_cell_ids - 1) {
                    printf(", ");
                }
            }
            printf("]\n");
        } else {
            printf("(not yet mapped)\n");
        }

        printf("\n");
    }

    free(indices);
}
