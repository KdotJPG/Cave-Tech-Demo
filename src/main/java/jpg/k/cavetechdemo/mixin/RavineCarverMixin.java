package jpg.k.cavetechdemo.mixin;

import net.minecraft.util.math.BlockPos;
import net.minecraft.world.biome.Biome;
import net.minecraft.world.chunk.Chunk;
import net.minecraft.world.gen.ProbabilityConfig;
import net.minecraft.world.gen.carver.RavineCarver;
import org.spongepowered.asm.mixin.Mixin;
import org.spongepowered.asm.mixin.Overwrite;

import java.util.BitSet;
import java.util.Random;
import java.util.function.Function;

@Mixin(RavineCarver.class)
public class RavineCarverMixin {

    /**
     * Disable Vanilla ravines, for this tech demo.
     * It would be better to either remove the carver from the biomes, or check here if it should run or not.
     * I will do that if this turns into any kind of released Minecraft mod.
     * @author K.jpg
     */
    @Overwrite
    public boolean carve(Chunk chunk, Function<BlockPos, Biome> posToBiome, Random random, int seaLevel, int j, int k, int l, int m, BitSet carvingMask, ProbabilityConfig probabilityConfig) {
        return false;
    }

}
