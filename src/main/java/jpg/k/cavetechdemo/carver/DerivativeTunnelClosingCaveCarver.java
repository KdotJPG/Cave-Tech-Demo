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
 * Derivative Tunnel Closing cave carver.
 *
 * This will calculate a sum of squares of two noise values to produce a roughly round tunnel shape.
 * Then, it will sample the derivative vectors and use them to close off parts that become too flat.
 * This avoids the least cave-like parts of this type of gen, and also introduces organic dead ends.
 *
 * Created by K.jpg on 5/27/2021.
 */
public class DerivativeTunnelClosingCaveCarver extends Carver<DefaultCarverConfig> {

    private static final double XZ_FREQUENCY_1 = 1.0 / 48.0;
    private static final double XZ_FREQUENCY_2 = 1.0 / 48.0;
    private static final double Y_FREQUENCY_1 = 1.0 / 32.0;
    private static final double Y_FREQUENCY_2 = 1.0 / 32.0;
    private static final double TUNNEL_THRESHOLD = 0.04;
    private static final double IDEAL_OPENSIMPLEX2_PAIR_VERTICAL_SAMPLING_OFFSET = 0.4330127018922193; // sqrt(3)/4, half the distance between two lattice points vertically
    private static final double CAVE_SHAPE_CLOSEOFF_SENSITIVITY = 2.0;

    // Use this simplex-type noise to encourage caves to follow a wide variety of directions.
    // This optimized version is also super fast for successive columnar evaluation.
    // In a real mod, we would need to pull in the world seed, instead of hardcoding it
    // 1.17 adds the world object to the chunk object, so it will be much easier there.
    OpenSimplex2CachedColumnNoise noise1 = new OpenSimplex2CachedColumnNoise(0);
    OpenSimplex2CachedColumnNoise noise2 = new OpenSimplex2CachedColumnNoise(1);

    public DerivativeTunnelClosingCaveCarver() {
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
        OpenSimplex2CachedColumnNoise.Values values1 = new OpenSimplex2CachedColumnNoise.Values();
        OpenSimplex2CachedColumnNoise.Values values2 = new OpenSimplex2CachedColumnNoise.Values();

        for (int z = 0; z < 16; z++) {
            int worldZ = (baseChunkZ << 4) | z;
            for (int x = 0; x < 16; x++) {
                int worldX = (baseChunkX << 4) | x;
                columnGen1.setXZ(worldX * XZ_FREQUENCY_1, worldZ * XZ_FREQUENCY_1);
                columnGen2.setXZ(worldX * XZ_FREQUENCY_2, worldZ * XZ_FREQUENCY_2);
                for (int y = heightmap.get(x, z); y >= 0; y--) {
                    double y1 = y * Y_FREQUENCY_1;
                    double y2 = y * Y_FREQUENCY_2 + IDEAL_OPENSIMPLEX2_PAIR_VERTICAL_SAMPLING_OFFSET;

                    // First, get one of the two noise values.
                    double value1 = columnGen1.getForY(y1);

                    // Start computing the sum of squares of what will be two noise values.
                    // This lets us often quickly rule out calculating the second noise.
                    double caveDensityValue = value1 * value1;
                    if (caveDensityValue >= TUNNEL_THRESHOLD) continue;

                    // Since it wasn't ruled out, compute the second noise value.
                    double value2 = columnGen2.getForY(y2);

                    // Finish computing the sum of their squares to create a tunnel shape that's roughly round.
                    caveDensityValue += value2 * value2;

                    // If this value doesn't pass the threshold test, we know there won't be a cave opening here.
                    if (caveDensityValue >= TUNNEL_THRESHOLD) continue;

                    // Now, we want to block out the parts of this cave gen where the carve-outs are more planar and less round.
                    // These occur when the noise values are flowing in almost the same direction. Get their derivative vectors.
                    // Getting only the value takes around 22ns on my laptop CPU, while this call takes around 34ns. This saves
                    // 12ns for most blocks in exchange for costing the extra 22ns for some blocks. It may be slightly better
                    // than that, too, considering these two calls never need to update the state of the noise column generator.
                    // It could also be optimized further by removing the value re-computation in this full evaluation call.
                    columnGen1.getForY(values1, y1);
                    columnGen2.getForY(values2, y2);

                    // Compute their squared magnitudes and dot product.
                    double noiseDeltaMagSq1 = values1.dx*values1.dx + values1.dy*values1.dy + values1.dz*values1.dz;
                    double noiseDeltaMagSq2 = values2.dx*values2.dx + values2.dy*values2.dy + values2.dz*values2.dz;
                    double noiseDeltaDot = values1.dx*values2.dx + values1.dy*values2.dy + values1.dz*values2.dz;

                    // If we normalized the derivative vectors (divided by the square roots of their magnitudes)
                    // then the dot product will be a value between -1 and 1. Closer to 0 would be when the noise
                    // values are flowing in roughly orthogonal directions (producing round caves), and closer to
                    // -1 or 1 would be when they're flowing in roughly parallel directions (producing flat parts).
                    // normalizedVector1 dot normalizedVector2 = (vector1 / |vector1|) dot (vector2 / |vector2|)
                    // = (vector1 dot vector2) / (|vector1| * |vector2|), which involves square roots.
                    // Square roots are good to avoid when possible, and there are often ways around them.
                    // Here, we want to square the dot product anyway, so that -1 and 1 both map to 1.
                    // Squaring the result makes the square roots in the denominator both disappear.
                    double squaredNormalizedDotProduct = (noiseDeltaDot * noiseDeltaDot) / (noiseDeltaMagSq1 * noiseDeltaMagSq2);

                    // As this gets closer to 1, we want to smoothly add to the cave density, so that it
                    // crosses the threshold before the flat areas come up. This will organically close off the
                    // caves while they're still round enough. Note that, since we're only adding to it, we aren't
                    // creating new areas that were skipped by the initial value check, but should be carved out now.
                    // Square the value again to make the close-off curve steeper, and prevent most isolated pockets
                    // that can appear following a close-off.
                    double caveClosingValue = (squaredNormalizedDotProduct * squaredNormalizedDotProduct) * (TUNNEL_THRESHOLD * CAVE_SHAPE_CLOSEOFF_SENSITIVITY);
                    caveDensityValue += caveClosingValue;

                    // Now test the updated density value against the threshold.
                    if (caveDensityValue >= TUNNEL_THRESHOLD) continue;

                    // This code produces the same result while avoiding the division. It didn't seem to be faster though.
                    // Perhaps this is because it needs three extra multiplications to replace only a single division.
                    /*double squaredNormalizedDotProductNumerator = noiseDeltaDot * noiseDeltaDot;
                    double squaredNormalizedDotProductDenominator = noiseDeltaMagSq1 * noiseDeltaMagSq2;
                    double twiceSquaredNormalizedDotProductDenominator = squaredNormalizedDotProductDenominator * squaredNormalizedDotProductDenominator;
                    double caveClosingValueNumerator = (squaredNormalizedDotProductNumerator * squaredNormalizedDotProductNumerator) * (TUNNEL_THRESHOLD * CAVE_SHAPE_CLOSEOFF_SENSITIVITY);
                    double caveDensityValueNumerator = caveDensityValue * twiceSquaredNormalizedDotProductDenominator + caveClosingValueNumerator;
                    if (caveDensityValueNumerator >= TUNNEL_THRESHOLD * twiceSquaredNormalizedDotProductDenominator) continue;*/

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
