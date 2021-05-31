package jpg.k.cavetechdemo.util.noise;

/**
 * OpenSimplex2 noise, speed-optimized for column generation.
 * Lock in an X/Z position, and call getForY(y) for each successive Y position.
 * It remembers enough of the current state of the noise to achieve a ~2.6x speedup.
 *
 * Works best for small steps in Y. Larger steps may see diminished performance gains.
 * Therefore it is more ideal a step towards avoiding interpolation, than to use with it.
 * Works great with conditional noise layer skipping!
 *
 * Uses "Improve XZ Planes" domain rotation to reduce subtle diagonal patterns.
 *
 * Not recommended for 2D noise where X and Y are both horizontal, as there is some bias.
 * I may later release noise better suited to this purpose.
 *
 * In this variant, I added support for derivatives.
 *
 * @author K.jpg
 */

public class OpenSimplex2CachedColumnNoise {

    private static final double ROOT3 = 1.7320508075688772;
    private static final double ROOT3OVER3 = 0.577350269189626;

    private static final int PERMUTATION_TABLE_SIZE_EXPONENT = 8;

    private static final int PSIZE = 1 << PERMUTATION_TABLE_SIZE_EXPONENT;
    private static final int PMASK = PSIZE - 1;

    private short[] perm;

    public OpenSimplex2CachedColumnNoise(long seed) {
        perm = new short[PSIZE];
        short[] source = new short[PSIZE];
        for (short i = 0; i < PSIZE; i++)
            source[i] = i;
        for (int i = PSIZE - 1; i >= 0; i--) {
            seed = seed * 6364136223846793005L + 1442695040888963407L;
            int r = (int)((seed + 31) % (i + 1));
            if (r < 0)
                r += (i + 1);
            perm[i] = source[r];
            source[r] = source[i];
        }
    }

    public ColumnGenerator columnGenerator() {
        return new ColumnGenerator();
    }

    public ColumnGeneratorWithDerivatives columnGeneratorWithDerivatives() {
        return new ColumnGeneratorWithDerivatives();
    }

    public class ColumnGenerator {

        private double xs, zs, xzr;
        private double localMinY, localMaxY;

        private double g0b, g1b, g2b, g3b;
        private double g0y, g1y, g2y, g3y;
        private double y0, y1, y2, y3;
        private double xzFalloff0, xzFalloff1, xzFalloff2, xzFalloff3;

        protected ColumnGenerator() {
            localMinY = Double.POSITIVE_INFINITY;
        }

        public void setXZ(double x, double z) {

            // Domain rotation, start
            double xz = x + z;
            double s2 = xz * -0.211324865405187;
            this.xs = x + s2;
            this.zs = z + s2;
            this.xzr = xz * -ROOT3OVER3;

            localMinY = Double.POSITIVE_INFINITY;
        }

        public double getForY(double y) {

            if (y < localMinY || y > localMaxY) {

                // Domain rotation, finish
                double yy = y * ROOT3OVER3;
                double xr = xs + yy;
                double zr = zs + yy;
                double yr = xzr + yy;

                // Grid base and bounds
                int xrb = fastRound(xr), yrb = fastRound(yr), zrb = fastRound(zr);
                double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

                // -1 if positive, 1 if negative
                int xNSign = (int)(-1.0 - xri) | 1;
                int yNSign = (int)(-1.0 - yri) | 1;
                int zNSign = (int)(-1.0 - zri) | 1;

                // Absolute values, using the above
                double axri = xNSign * -xri;
                double ayri = yNSign * -yri;
                double azri = zNSign * -zri;

                for (int l = 0; ; l++)
                {
                    // Get delta x and z in world-space coordinates away from grid base.
                    // This way we can avoid re-computing the stuff for X and Z, and
                    // only update for Y as Y updates.
                    double s2i = (xri + zri) * -0.211324865405187 - (yri * ROOT3OVER3);
                    double xi0 = xri + s2i;
                    double zi0 = zri + s2i;
                    double yi0 = (xri + zri) * ROOT3OVER3 + (yri * ROOT3OVER3);

                    // Closest vertex on this half-grid
                    int gi = gradIndex(xrb, yrb, zrb);

                    // Its gradient and falloff data
                    g0b = RGRADIENTS_3D[gi | 0] * xi0 + RGRADIENTS_3D[gi | 1] * zi0; g0y = RGRADIENTS_3D[gi | 2];
                    xzFalloff0 = 0.6 - xi0 * xi0 - zi0 * zi0;
                    y0 = y - yi0;

                    // Find second-closest vertex
                    if (axri >= ayri && axri >= azri) {
                        gi = gradIndex(xrb - xNSign, yrb, zrb);

                        // Position of this vertex, in world-space coordinates
                        double xi = xi0 + xNSign * 0.788675134594813, zi = zi0 + xNSign * -0.211324865405187, yi = yi0 + xNSign * ROOT3OVER3;

                        // Gradient and falloff data
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;

                        // One of the planar boundaries where the state of the noise changes lies between the two
                        // diagonally-opposite vertices on the other half-grid which connect to that half-grid's
                        // closest vertex via edges perpendicular to the edge formed on this half-grid.
                        // We store to localMinY, but we don't yet know if it's the min or the max.
                        // We'll resolve that later.
                        localMinY = y - (ROOT3 / 2) * (yri + zri);
                    } else if (ayri > axri && ayri >= azri) {
                        gi = gradIndex(xrb, yrb - yNSign, zrb);
                        double xi = xi0 + yNSign * -ROOT3OVER3, zi = zi0 + yNSign * -ROOT3OVER3, yi = yi0 + yNSign * ROOT3OVER3;
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;
                        localMinY = y - (ROOT3 / 2) * (xri + zri);
                    } else {
                        gi = gradIndex(xrb, yrb, zrb - zNSign);
                        double xi = xi0 + zNSign * -0.211324865405187, zi = zi0 + zNSign * 0.788675134594813, yi = yi0 + zNSign * ROOT3OVER3;
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;
                        localMinY = y - (ROOT3 / 2) * (xri + yri);
                    }

                    if (l == 1) break;

                    // Flip everyhing to reference the closest vertex on the other half-grid
                    axri = 0.5f - axri;
                    ayri = 0.5f - ayri;
                    azri = 0.5f - azri;
                    xri = xNSign * axri;
                    yri = yNSign * ayri;
                    zri = zNSign * azri;
                    xrb -= (xNSign >> 1) - (PSIZE / 2);
                    yrb -= (yNSign >> 1) - (PSIZE / 2);
                    zrb -= (zNSign >> 1) - (PSIZE / 2);
                    xNSign = -xNSign;
                    yNSign = -yNSign;
                    zNSign = -zNSign;

                    // Shift these over to vacate them for the other iteration.
                    g2b = g0b;
                    g3b = g1b;
                    g2y = g0y;
                    g3y = g1y;
                    y2 = y0;
                    y3 = y1;
                    xzFalloff2 = xzFalloff0;
                    xzFalloff3 = xzFalloff1;
                    localMaxY = localMinY;
                }

                // It isn't guaranteed which will be the min or max, so correct that here.
                if (localMinY > localMaxY) {
                    double temp = localMaxY;
                    localMaxY = localMinY;
                    localMinY = temp;
                }
            }
            double value = 0;

            double dy0 = y - y0;
            double dy0Sq = dy0 * dy0;
            if (xzFalloff0 > dy0Sq) {
                double falloff = xzFalloff0 - dy0Sq;
                falloff *= falloff;
                value = falloff * falloff * (g0b + g0y * dy0);
            }

            double dy1 = y - y1;
            double dy1Sq = dy1 * dy1;
            if (xzFalloff1 > dy1Sq) {
                double falloff = xzFalloff1 - dy1Sq;
                falloff *= falloff;
                value += falloff * falloff * (g1b + g1y * dy1);
            }

            double dy2 = y - y2;
            double dy2Sq = dy2 * dy2;
            if (xzFalloff2 > dy2Sq) {
                double falloff = xzFalloff2 - dy2Sq;
                falloff *= falloff;
                value += falloff * falloff * (g2b + g2y * dy2);
            }

            double dy3 = y - y3;
            double dy3Sq = dy3 * dy3;
            if (xzFalloff3 > dy3Sq) {
                double falloff = xzFalloff3 - dy3Sq;
                falloff *= falloff;
                value += falloff * falloff * (g3b + g3y * dy3);
            }

            return value;
        }

    }

    public class ColumnGeneratorWithDerivatives {

        private double xs, zs, xzr;
        private double localMinY, localMaxY;

        private double g0b, g0x, g0y, g0z;
        private double g1b, g1x, g1y, g1z;
        private double g2b, g2x, g2y, g2z;
        private double g3b, g3x, g3y, g3z;
        private double dx0, dz0, dx1, dz1;
        private double dx2, dz2, dx3, dz3;

        private double y0, y1, y2, y3;
        private double xzFalloff0, xzFalloff1, xzFalloff2, xzFalloff3;

        protected ColumnGeneratorWithDerivatives() {
            localMinY = Double.POSITIVE_INFINITY;
        }

        public void setXZ(double x, double z) {

            // Domain rotation, start
            double xz = x + z;
            double s2 = xz * -0.211324865405187;
            this.xs = x + s2;
            this.zs = z + s2;
            this.xzr = xz * -ROOT3OVER3;

            localMinY = Double.POSITIVE_INFINITY;
        }

        public void getForY(Values values, double y) {

            if (y < localMinY || y > localMaxY) {

                // Domain rotation, finish
                double yy = y * ROOT3OVER3;
                double xr = xs + yy;
                double zr = zs + yy;
                double yr = xzr + yy;

                // Grid base and bounds
                int xrb = fastRound(xr), yrb = fastRound(yr), zrb = fastRound(zr);
                double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

                // -1 if positive, 1 if negative
                int xNSign = (int)(-1.0 - xri) | 1;
                int yNSign = (int)(-1.0 - yri) | 1;
                int zNSign = (int)(-1.0 - zri) | 1;

                // Absolute values, using the above
                double axri = xNSign * -xri;
                double ayri = yNSign * -yri;
                double azri = zNSign * -zri;

                for (int l = 0; ; l++)
                {
                    // Get delta x and z in world-space coordinates away from grid base.
                    // This way we can avoid re-computing the stuff for X and Z, and
                    // only update for Y as Y updates.
                    double s2i = (xri + zri) * -0.211324865405187 - (yri * ROOT3OVER3);
                    double xi0 = xri + s2i;
                    double zi0 = zri + s2i;
                    double yi0 = (xri + zri) * ROOT3OVER3 + (yri * ROOT3OVER3);
                    this.dx0 = xi0; this.dz0 = zi0;

                    // Closest vertex on this half-grid
                    int gi = gradIndex(xrb, yrb, zrb);

                    // Its gradient and falloff data
                    g0b = RGRADIENTS_3D[gi | 0] * xi0 + RGRADIENTS_3D[gi | 1] * zi0;
                    g0x = RGRADIENTS_3D[gi | 0]; g0z = RGRADIENTS_3D[gi | 1]; g0y = RGRADIENTS_3D[gi | 2];
                    xzFalloff0 = 0.6 - xi0 * xi0 - zi0 * zi0;
                    y0 = y - yi0;

                    // Find second-closest vertex
                    if (axri >= ayri && axri >= azri) {
                        gi = gradIndex(xrb - xNSign, yrb, zrb);

                        // Position of this vertex, in world-space coordinates
                        double xi = xi0 + xNSign * 0.788675134594813, zi = zi0 + xNSign * -0.211324865405187, yi = yi0 + xNSign * ROOT3OVER3;
                        this.dx1 = xi; this.dz1 = zi;

                        // Gradient and falloff data
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                        g1x = RGRADIENTS_3D[gi | 0]; g1z = RGRADIENTS_3D[gi | 1]; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;

                        // One of the planar boundaries where the state of the noise changes lies between the two
                        // diagonally-opposite vertices on the other half-grid which connect to that half-grid's
                        // closest vertex via edges perpendicular to the edge formed on this half-grid.
                        // We store to localMinY, but we don't yet know if it's the min or the max.
                        // We'll resolve that later.
                        localMinY = y - (ROOT3 / 2) * (yri + zri);
                    } else if (ayri > axri && ayri >= azri) {
                        gi = gradIndex(xrb, yrb - yNSign, zrb);
                        double xi = xi0 + yNSign * -ROOT3OVER3, zi = zi0 + yNSign * -ROOT3OVER3, yi = yi0 + yNSign * ROOT3OVER3;
                        this.dx1 = xi; this.dz1 = zi;
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                        g1x = RGRADIENTS_3D[gi | 0]; g1z = RGRADIENTS_3D[gi | 1]; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;
                        localMinY = y - (ROOT3 / 2) * (xri + zri);
                    } else {
                        gi = gradIndex(xrb, yrb, zrb - zNSign);
                        double xi = xi0 + zNSign * -0.211324865405187, zi = zi0 + zNSign * 0.788675134594813, yi = yi0 + zNSign * ROOT3OVER3;
                        this.dx1 = xi; this.dz1 = zi;
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                        g1x = RGRADIENTS_3D[gi | 0]; g1z = RGRADIENTS_3D[gi | 1]; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;
                        localMinY = y - (ROOT3 / 2) * (xri + yri);
                    }

                    if (l == 1) break;

                    // Flip everyhing to reference the closest vertex on the other half-grid
                    axri = 0.5f - axri;
                    ayri = 0.5f - ayri;
                    azri = 0.5f - azri;
                    xri = xNSign * axri;
                    yri = yNSign * ayri;
                    zri = zNSign * azri;
                    xrb -= (xNSign >> 1) - (PSIZE / 2);
                    yrb -= (yNSign >> 1) - (PSIZE / 2);
                    zrb -= (zNSign >> 1) - (PSIZE / 2);
                    xNSign = -xNSign;
                    yNSign = -yNSign;
                    zNSign = -zNSign;

                    // Shift these over to vacate them for the other iteration.
                    g2b = g0b;
                    g3b = g1b;
                    g2x = g0x;
                    g3x = g1x;
                    g2z = g0z;
                    g3z = g1z;
                    g2y = g0y;
                    g3y = g1y;
                    dx2 = dx0;
                    dx3 = dx1;
                    dz2 = dz0;
                    dz3 = dz1;
                    y2 = y0;
                    y3 = y1;
                    xzFalloff2 = xzFalloff0;
                    xzFalloff3 = xzFalloff1;
                    localMaxY = localMinY;
                }

                // It isn't guaranteed which will be the min or max, so correct that here.
                if (localMinY > localMaxY) {
                    double temp = localMaxY;
                    localMaxY = localMinY;
                    localMinY = temp;
                }
            }
            double value = 0, dx = 0, dy = 0, dz = 0;

            double dy0 = y - y0;
            double dy0Sq = dy0 * dy0;
            if (xzFalloff0 > dy0Sq) {
                double falloff = xzFalloff0 - dy0Sq;
                double falloffSq = falloff * falloff;
                double ramp = g0b + g0y * dy0;
                value = falloffSq * falloffSq * ramp;
                dx = (g0x * falloff - 8 * ramp * dx0) * (falloff * falloffSq);
                dy = (g0y * falloff - 8 * ramp * dy0) * (falloff * falloffSq);
                dz = (g0z * falloff - 8 * ramp * dz0) * (falloff * falloffSq);
            }

            double dy1 = y - y1;
            double dy1Sq = dy1 * dy1;
            if (xzFalloff1 > dy1Sq) {
                double falloff = xzFalloff1 - dy1Sq;
                double falloffSq = falloff * falloff;
                double ramp = g1b + g1y * dy1;
                value += falloffSq * falloffSq * ramp;
                dx += (g1x * falloff - 8 * ramp * dx1) * (falloff * falloffSq);
                dy += (g1y * falloff - 8 * ramp * dy1) * (falloff * falloffSq);
                dz += (g1z * falloff - 8 * ramp * dz1) * (falloff * falloffSq);
            }

            double dy2 = y - y2;
            double dy2Sq = dy2 * dy2;
            if (xzFalloff2 > dy2Sq) {
                double falloff = xzFalloff2 - dy2Sq;
                double falloffSq = falloff * falloff;
                double ramp = g2b + g2y * dy2;
                value += falloffSq * falloffSq * ramp;
                dx += (g2x * falloff - 8 * ramp * dx2) * (falloff * falloffSq);
                dy += (g2y * falloff - 8 * ramp * dy2) * (falloff * falloffSq);
                dz += (g2z * falloff - 8 * ramp * dz2) * (falloff * falloffSq);
            }

            double dy3 = y - y3;
            double dy3Sq = dy3 * dy3;
            if (xzFalloff3 > dy3Sq) {
                double falloff = xzFalloff3 - dy3Sq;
                double falloffSq = falloff * falloff;
                double ramp = g3b + g3y * dy3;
                value += falloffSq * falloffSq * ramp;
                dx += (g3x * falloff - 8 * ramp * dx3) * (falloff * falloffSq);
                dy += (g3y * falloff - 8 * ramp * dy3) * (falloff * falloffSq);
                dz += (g3z * falloff - 8 * ramp * dz3) * (falloff * falloffSq);
            }

            values.value = value;
            values.dx = dx;
            values.dy = dy;
            values.dz = dz;
        }

        public double getForY(double y) {

            if (y < localMinY || y > localMaxY) {

                // Domain rotation, finish
                double yy = y * ROOT3OVER3;
                double xr = xs + yy;
                double zr = zs + yy;
                double yr = xzr + yy;

                // Grid base and bounds
                int xrb = fastRound(xr), yrb = fastRound(yr), zrb = fastRound(zr);
                double xri = xr - xrb, yri = yr - yrb, zri = zr - zrb;

                // -1 if positive, 1 if negative
                int xNSign = (int)(-1.0 - xri) | 1;
                int yNSign = (int)(-1.0 - yri) | 1;
                int zNSign = (int)(-1.0 - zri) | 1;

                // Absolute values, using the above
                double axri = xNSign * -xri;
                double ayri = yNSign * -yri;
                double azri = zNSign * -zri;

                for (int l = 0; ; l++)
                {
                    // Get delta x and z in world-space coordinates away from grid base.
                    // This way we can avoid re-computing the stuff for X and Z, and
                    // only update for Y as Y updates.
                    double s2i = (xri + zri) * -0.211324865405187 - (yri * ROOT3OVER3);
                    double xi0 = xri + s2i;
                    double zi0 = zri + s2i;
                    double yi0 = (xri + zri) * ROOT3OVER3 + (yri * ROOT3OVER3);
                    this.dx0 = xi0; this.dz0 = zi0;

                    // Closest vertex on this half-grid
                    int gi = gradIndex(xrb, yrb, zrb);

                    // Its gradient and falloff data
                    g0b = RGRADIENTS_3D[gi | 0] * xi0 + RGRADIENTS_3D[gi | 1] * zi0;
                    g0x = RGRADIENTS_3D[gi | 0]; g0z = RGRADIENTS_3D[gi | 1]; g0y = RGRADIENTS_3D[gi | 2];
                    xzFalloff0 = 0.6 - xi0 * xi0 - zi0 * zi0;
                    y0 = y - yi0;

                    // Find second-closest vertex
                    if (axri >= ayri && axri >= azri) {
                        gi = gradIndex(xrb - xNSign, yrb, zrb);

                        // Position of this vertex, in world-space coordinates
                        double xi = xi0 + xNSign * 0.788675134594813, zi = zi0 + xNSign * -0.211324865405187, yi = yi0 + xNSign * ROOT3OVER3;
                        this.dx1 = xi; this.dz1 = zi;

                        // Gradient and falloff data
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                        g1x = RGRADIENTS_3D[gi | 0]; g1z = RGRADIENTS_3D[gi | 1]; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;

                        // One of the planar boundaries where the state of the noise changes lies between the two
                        // diagonally-opposite vertices on the other half-grid which connect to that half-grid's
                        // closest vertex via edges perpendicular to the edge formed on this half-grid.
                        // We store to localMinY, but we don't yet know if it's the min or the max.
                        // We'll resolve that later.
                        localMinY = y - (ROOT3 / 2) * (yri + zri);
                    } else if (ayri > axri && ayri >= azri) {
                        gi = gradIndex(xrb, yrb - yNSign, zrb);
                        double xi = xi0 + yNSign * -ROOT3OVER3, zi = zi0 + yNSign * -ROOT3OVER3, yi = yi0 + yNSign * ROOT3OVER3;
                        this.dx1 = xi; this.dz1 = zi;
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                        g1x = RGRADIENTS_3D[gi | 0]; g1z = RGRADIENTS_3D[gi | 1]; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;
                        localMinY = y - (ROOT3 / 2) * (xri + zri);
                    } else {
                        gi = gradIndex(xrb, yrb, zrb - zNSign);
                        double xi = xi0 + zNSign * -0.211324865405187, zi = zi0 + zNSign * 0.788675134594813, yi = yi0 + zNSign * ROOT3OVER3;
                        this.dx1 = xi; this.dz1 = zi;
                        g1b = RGRADIENTS_3D[gi | 0] * xi + RGRADIENTS_3D[gi | 1] * zi;
                        g1x = RGRADIENTS_3D[gi | 0]; g1z = RGRADIENTS_3D[gi | 1]; g1y = RGRADIENTS_3D[gi | 2];
                        xzFalloff1 = 0.6 - xi * xi - zi * zi;
                        y1 = y - yi;
                        localMinY = y - (ROOT3 / 2) * (xri + yri);
                    }

                    if (l == 1) break;

                    // Flip everyhing to reference the closest vertex on the other half-grid
                    axri = 0.5f - axri;
                    ayri = 0.5f - ayri;
                    azri = 0.5f - azri;
                    xri = xNSign * axri;
                    yri = yNSign * ayri;
                    zri = zNSign * azri;
                    xrb -= (xNSign >> 1) - (PSIZE / 2);
                    yrb -= (yNSign >> 1) - (PSIZE / 2);
                    zrb -= (zNSign >> 1) - (PSIZE / 2);
                    xNSign = -xNSign;
                    yNSign = -yNSign;
                    zNSign = -zNSign;

                    // Shift these over to vacate them for the other iteration.
                    g2b = g0b;
                    g3b = g1b;
                    g2x = g0x;
                    g3x = g1x;
                    g2z = g0z;
                    g3z = g1z;
                    g2y = g0y;
                    g3y = g1y;
                    dx2 = dx0;
                    dx3 = dx1;
                    dz2 = dz0;
                    dz3 = dz1;
                    y2 = y0;
                    y3 = y1;
                    xzFalloff2 = xzFalloff0;
                    xzFalloff3 = xzFalloff1;
                    localMaxY = localMinY;
                }

                // It isn't guaranteed which will be the min or max, so correct that here.
                if (localMinY > localMaxY) {
                    double temp = localMaxY;
                    localMaxY = localMinY;
                    localMinY = temp;
                }
            }
            double value = 0, dx = 0, dy = 0, dz = 0;

            double dy0 = y - y0;
            double dy0Sq = dy0 * dy0;
            if (xzFalloff0 > dy0Sq) {
                double falloff = xzFalloff0 - dy0Sq;
                double falloffSq = falloff * falloff;
                double ramp = g0b + g0y * dy0;
                value = falloffSq * falloffSq * ramp;
            }

            double dy1 = y - y1;
            double dy1Sq = dy1 * dy1;
            if (xzFalloff1 > dy1Sq) {
                double falloff = xzFalloff1 - dy1Sq;
                double falloffSq = falloff * falloff;
                double ramp = g1b + g1y * dy1;
                value += falloffSq * falloffSq * ramp;
            }

            double dy2 = y - y2;
            double dy2Sq = dy2 * dy2;
            if (xzFalloff2 > dy2Sq) {
                double falloff = xzFalloff2 - dy2Sq;
                double falloffSq = falloff * falloff;
                double ramp = g2b + g2y * dy2;
                value += falloffSq * falloffSq * ramp;
            }

            double dy3 = y - y3;
            double dy3Sq = dy3 * dy3;
            if (xzFalloff3 > dy3Sq) {
                double falloff = xzFalloff3 - dy3Sq;
                double falloffSq = falloff * falloff;
                double ramp = g3b + g3y * dy3;
                value += falloffSq * falloffSq * ramp;
            }

            return value;
        }

    }

    private static int fastRound(double x) {
        return (int)(x >= 0 ? x + 0.5 : x - 0.5);
    }

    private int gradIndex(int xrv, int yrv, int zrv) {
        return perm[perm[perm[xrv & PMASK] ^ (yrv & PMASK)] ^ (zrv & PMASK)] & 0xFC;
    }

    private static final double[] RGRADIENTS_3D = {
            0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
            1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
            1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
            0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
            1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
            1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
            0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
            1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
            1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
            0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
            1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
            1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
            0, 1, 1, 0,  0,-1, 1, 0,  0, 1,-1, 0,  0,-1,-1, 0,
            1, 0, 1, 0, -1, 0, 1, 0,  1, 0,-1, 0, -1, 0,-1, 0,
            1, 1, 0, 0, -1, 1, 0, 0,  1,-1, 0, 0, -1,-1, 0, 0,
            1, 1, 0, 0,  0,-1, 1, 0, -1, 1, 0, 0,  0,-1,-1, 0
    };

    private static final double NORMALIZING_MULTIPLIER_3D = 32.69428253173828125;
    static {
        for (int i = 0; i < RGRADIENTS_3D.length; i += 4) {
            double gx = RGRADIENTS_3D[i | 0];
            double gy = RGRADIENTS_3D[i | 1];
            double gz = RGRADIENTS_3D[i | 2];

            // Rotate and apply noise normalizing constant
            // Store as XZY so X and Z are next to each other when we read them at the same time.
            double s2 = (gx + gz) * -0.211324865405187 + gy * 0.577350269189626;
            RGRADIENTS_3D[i | 0] = (gx + s2) * NORMALIZING_MULTIPLIER_3D;
            RGRADIENTS_3D[i | 1] = (gz + s2) * NORMALIZING_MULTIPLIER_3D;
            RGRADIENTS_3D[i | 2] = (gy - gx - gz) * (ROOT3OVER3 * NORMALIZING_MULTIPLIER_3D);
        }
    }

    public static class Values {
        public double value, dx, dy, dz;
    }

}