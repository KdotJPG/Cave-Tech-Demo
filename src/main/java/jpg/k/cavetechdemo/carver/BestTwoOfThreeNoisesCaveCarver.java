package jpg.k.cavetechdemo.carver;

import jpg.k.cavetechdemo.util.noise.OpenSimplex2CachedColumnNoise;
import net.minecraft.util.math.BlockPos;
import net.minecraft.world.Heightmap;
import net.minecraft.world.biome.Biome;
import net.minecraft.world.chunk.Chunk;
import net.minecraft.world.gen.carver.Carver;
import net.minecraft.world.gen.carver.DefaultCarverConfig;
import org.apache.commons.lang3.mutable.MutableBoolean;

import java.util.BitSet;
import java.util.Random;
import java.util.function.Function;

/**
 * Best-Two-Of-Three-Noises cave carver.
 *
 * This was a bit of a failed experiment. My goal was to smoothly and dynamically pick the best two of three noises
 * to carve caves with, but it didn't end up producing caves the way I intended. It's possible I made a math error,
 * but it's also possible the idea just doesn't work out as I had hoped. Leaving this in here in case anyone finds
 * parts of it useful. Particularly, the smooth weight selection formula might be great for other effects.
 *
 * Created by K.jpg on 5/31/2021.
 */
public class BestTwoOfThreeNoisesCaveCarver extends Carver<DefaultCarverConfig> {

    private static final double XZ_FREQUENCY_1 = 1.0 / 48.0;
    private static final double XZ_FREQUENCY_2 = 1.0 / 48.0;
    private static final double XZ_FREQUENCY_3 = 1.0 / 48.0;
    private static final double Y_FREQUENCY_1 = 1.0 / 32.0;
    private static final double Y_FREQUENCY_2 = 1.0 / 32.0;
    private static final double Y_FREQUENCY_3 = 1.0 / 32.0;
    private static final double TUNNEL_THRESHOLD = 0.04;
    private static final double TUNNEL_THRESHOLD_SQRT = Math.sqrt(TUNNEL_THRESHOLD);
    private static final double IDEAL_OPENSIMPLEX2_TRIPLET_VERTICAL_SAMPLING_OFFSET = 0.28867513459481287; // sqrt(3)/6, one third the distance between two lattice points vertically
    private static final double CAVE_NOISE_TRANSITION = 0.2;

    // Use this simplex-type noise to encourage caves to follow a wide variety of directions.
    // This optimized version is also super fast for successive columnar evaluation.
    // In a real mod, we would need to pull in the world seed, instead of hardcoding it
    // 1.17 adds the world object to the chunk object, so it will be much easier there.
    OpenSimplex2CachedColumnNoise noise1 = new OpenSimplex2CachedColumnNoise(0);
    OpenSimplex2CachedColumnNoise noise2 = new OpenSimplex2CachedColumnNoise(1);
    OpenSimplex2CachedColumnNoise noise3 = new OpenSimplex2CachedColumnNoise(2);

    public BestTwoOfThreeNoisesCaveCarver() {
        super(DefaultCarverConfig.CODEC, 256);
    }

    @Override
    public boolean carve(Chunk chunk, Function<BlockPos, Biome> posToBiome, Random random, int seaLevel, int targetChunkX, int targetChunkZ, int baseChunkX, int baseChunkZ, BitSet carvingMask, DefaultCarverConfig config) {
        if (targetChunkX != baseChunkX || targetChunkZ != baseChunkZ) return true; // Should I be returning true or false here?

        long start = System.currentTimeMillis();

        BlockPos.Mutable mutable = new BlockPos.Mutable();
        BlockPos.Mutable mutable2 = new BlockPos.Mutable();
        BlockPos.Mutable mutable3 = new BlockPos.Mutable();
        MutableBoolean mutableBoolean = new MutableBoolean(false);
        Heightmap heightmap = chunk.getHeightmap(Heightmap.Type.WORLD_SURFACE_WG);

        OpenSimplex2CachedColumnNoise.ColumnGeneratorWithDerivatives columnGen1 = noise1.columnGeneratorWithDerivatives();
        OpenSimplex2CachedColumnNoise.ColumnGeneratorWithDerivatives columnGen2 = noise2.columnGeneratorWithDerivatives();
        OpenSimplex2CachedColumnNoise.ColumnGeneratorWithDerivatives columnGen3 = noise3.columnGeneratorWithDerivatives();
        OpenSimplex2CachedColumnNoise.Values values1 = new OpenSimplex2CachedColumnNoise.Values();
        OpenSimplex2CachedColumnNoise.Values values2 = new OpenSimplex2CachedColumnNoise.Values();
        OpenSimplex2CachedColumnNoise.Values values3 = new OpenSimplex2CachedColumnNoise.Values();

        for (int z = 0; z < 16; z++) {
            int worldZ = (baseChunkZ << 4) | z;
            for (int x = 0; x < 16; x++) {
                int worldX = (baseChunkX << 4) | x;
                columnGen1.setXZ(worldX * XZ_FREQUENCY_1, worldZ * XZ_FREQUENCY_1);
                columnGen2.setXZ(worldX * XZ_FREQUENCY_2, worldZ * XZ_FREQUENCY_2);
                columnGen3.setXZ(worldX * XZ_FREQUENCY_3, worldZ * XZ_FREQUENCY_3);
                for (int y = heightmap.get(x, z); y >= 0; y--) {
                    double y1 = y * Y_FREQUENCY_1;
                    double y2 = y * Y_FREQUENCY_2 + IDEAL_OPENSIMPLEX2_TRIPLET_VERTICAL_SAMPLING_OFFSET;
                    double y3 = y * Y_FREQUENCY_2 + IDEAL_OPENSIMPLEX2_TRIPLET_VERTICAL_SAMPLING_OFFSET * 2;

                    // Get first two noise values.
                    double value1 = columnGen1.getForY(y1);
                    double value2 = columnGen2.getForY(y2);

                    // Given that at least one of these will be involved in the final density formula,
                    // it is impossible for there to be a cave here if both of them already exceed the range.
                    if ((value1 < -TUNNEL_THRESHOLD_SQRT || value1 > TUNNEL_THRESHOLD_SQRT)
                            && (value2 < -TUNNEL_THRESHOLD_SQRT || value2 > TUNNEL_THRESHOLD_SQRT))
                        continue;

                    // Get the third value.
                    double value3 = columnGen2.getForY(y3);

                    // The same can be said for third value combined with any of the previous.
                    if (value3 < -TUNNEL_THRESHOLD_SQRT || value3 > TUNNEL_THRESHOLD_SQRT) {
                        if (value1 < -TUNNEL_THRESHOLD_SQRT || value1 > TUNNEL_THRESHOLD_SQRT) continue;
                        if (value2 < -TUNNEL_THRESHOLD_SQRT || value2 > TUNNEL_THRESHOLD_SQRT) continue;
                    }

                    // Now get the derivatives (more expensive operation).
                    columnGen1.getForY(values1, y1);
                    columnGen2.getForY(values2, y2);
                    columnGen3.getForY(values3, y3);

                    // Compute squared normalized dot products of derivatives of pairs of noises.
                    double noiseDeltaMagSq1 = values1.dx*values1.dx + values1.dy*values1.dy + values1.dz*values1.dz;
                    double noiseDeltaMagSq2 = values2.dx*values2.dx + values2.dy*values2.dy + values2.dz*values2.dz;
                    double noiseDeltaMagSq3 = values3.dx*values3.dx + values3.dy*values3.dy + values3.dz*values3.dz;
                    double noiseDeltaDot12 = values1.dx*values2.dx + values1.dy*values2.dy + values1.dz*values2.dz;
                    double noiseDeltaDot13 = values1.dx*values3.dx + values1.dy*values3.dy + values1.dz*values3.dz;
                    double noiseDeltaDot23 = values2.dx*values3.dx + values2.dy*values3.dy + values2.dz*values3.dz;
                    double squaredNormalizedDot12 = (noiseDeltaDot12 * noiseDeltaDot12) / (noiseDeltaMagSq1 * noiseDeltaMagSq2);
                    double squaredNormalizedDot13 = (noiseDeltaDot13 * noiseDeltaDot13) / (noiseDeltaMagSq1 * noiseDeltaMagSq3);
                    double squaredNormalizedDot23 = (noiseDeltaDot23 * noiseDeltaDot23) / (noiseDeltaMagSq2 * noiseDeltaMagSq3);

                    // Fade curves between pairs of normalized squared dot products within the transition zones.
                    // These will be used to fade between which pair of noises is used based on how close their
                    // dot product is to zero.
                    double fade12 = smoothedClampedFadeWideDomain((squaredNormalizedDot23 - squaredNormalizedDot13) * (2.0 / CAVE_NOISE_TRANSITION));
                    double fade13 = smoothedClampedFadeWideDomain((squaredNormalizedDot23 - squaredNormalizedDot12) * (2.0 / CAVE_NOISE_TRANSITION));
                    double fade23 = smoothedClampedFadeWideDomain((squaredNormalizedDot13 - squaredNormalizedDot12) * (2.0 / CAVE_NOISE_TRANSITION));

                    // These weights (smoothly) pick which noise has the highest dot product with the other two,
                    // to (smoothly) remove it from the formula and let the other two form the caves.
                    double weight1For23 = (1 - fade12) * (1 - fade13);
                    double weight2For13 = fade12 * (1 - fade23);
                    double weight3For12 = fade23 * fade13;

                    // Normalize weights.
                    // Most of the time, the weight sum will be exactly 1. This is because the clamped fade curves will
                    // output 1 or 0 a lot, and they will be multiplied together. Thus it is meaningful to do an exact
                    // comparison on a floating point type, which is not ordinarily done. Checking this, we can skip.
                    double weightSum = weight1For23 + weight2For13 + weight3For12;
                    if (weightSum != 1.0) {
                        double inverseWeightSum = 1.0 / weightSum;
                        weight1For23 *= inverseWeightSum;
                        weight2For13 *= inverseWeightSum;
                        weight3For12 *= inverseWeightSum;
                    }

                    //double density = weight1For23 * (value2 * value2 + value3 * value3) + weight2For13 * (value1 * value1 + value3 * value3) + weight3For12 * (value1 * value1 + value2 * value2);
                    double density = (value1 * value1) * (weight2For13 + weight3For12) + (value2 * value2) * (weight1For23 + weight3For12) + (value3 * value3) * (weight1For23 + weight2For13);

                    // Now test the fully computed density value against the threshold.
                    if (density >= TUNNEL_THRESHOLD) continue;

                    // If all is good, we can carve out the block.
                    super.carveAtPoint(chunk, posToBiome, carvingMask, random, mutable, mutable2, mutable3, seaLevel, baseChunkX, baseChunkZ, worldX, worldZ, x, y, z, mutableBoolean);
                }
            }
        }

        long elapsed = System.currentTimeMillis() - start;
        recordChunkRuntime(elapsed);

        return true;
    }

    private long totalElapsed = 0;
    private int totalNumChunks = 0;
    private synchronized void recordChunkRuntime(long elapsed) {
        totalElapsed += elapsed;
        totalNumChunks++;
        if ((totalNumChunks & 0x1F) == 0) {
            System.out.println("Cave Generation milliseconds per chunk: " + totalElapsed * (1.0 / totalNumChunks));
            //totalElapsed = 0;
            //totalNumChunks = 0;
        }
    }

    // Input range clamped to -1 to 1
    // then mapped into the range 0 to 1
    // smoothed at the boundaries
    // https://www.desmos.com/calculator/yfshhrqmmg
    private double smoothedClampedFadeWideDomain(double t) {
        if (t <= -1) return 0;
        else if (t >= 1) return 1;
        else return 0.5 + t * (0.75 - t * t * 0.25);
    }

    @Override
    public boolean shouldCarve(Random var1, int var2, int var3, DefaultCarverConfig var4) {
        return true;
    }

    @Override
    protected boolean isPositionExcluded(double var1, double var3, double var5, int var7)
    {
        return true;
    }


}
