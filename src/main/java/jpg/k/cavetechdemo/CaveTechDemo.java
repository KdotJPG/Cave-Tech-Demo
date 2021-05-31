package jpg.k.cavetechdemo;

import net.fabricmc.api.ModInitializer;
import net.fabricmc.fabric.api.biome.v1.BiomeModifications;
import net.fabricmc.fabric.api.biome.v1.BiomeSelectors;
import net.minecraft.util.Identifier;
import net.minecraft.util.registry.Registry;

import jpg.k.cavetechdemo.carver.*;
import net.minecraft.util.registry.RegistryKey;
import net.minecraft.util.registry.BuiltinRegistries;
import net.minecraft.world.gen.GenerationStep;
import net.minecraft.world.gen.carver.Carver;
import net.minecraft.world.gen.carver.CarverConfig;
import net.minecraft.world.gen.carver.ConfiguredCarver;
import net.minecraft.world.gen.carver.DefaultCarverConfig;

public class CaveTechDemo implements ModInitializer {

    private static final DerivativeTunnelClosingCaveCarver DERIVATIVE_TUNNEL_CLOSING_CAVE_CARVER = new DerivativeTunnelClosingCaveCarver();
    private static final BestTwoOfThreeNoisesCaveCarver BEST_TWO_OF_THREE_NOISES_CAVE_CARVER = new BestTwoOfThreeNoisesCaveCarver();

	@Override
	public void onInitialize() {

	    // Change this constant to change which carver gets used in this tech demo.
        Identifier carverId = registerCarver(DERIVATIVE_TUNNEL_CLOSING_CAVE_CARVER, "demo_carver");

		// BiomeModifications is an experimental Fabric API feature. Works fine for this demo.
        BiomeModifications.addCarver(BiomeSelectors.all(), GenerationStep.Carver.AIR, RegistryKey.of(Registry.CONFIGURED_CARVER_WORLDGEN, carverId));

    }

    private <T extends CarverConfig> Identifier registerCarver(Carver<T> carver, String carverIdPathName) {
        Identifier carverId = new Identifier("cavetechdemo", "demo_carver");
        ConfiguredCarver<T> configuredCarver = new ConfiguredCarver<>(carver, (T)CarverConfig.DEFAULT);
        Registry.register(Registry.CARVER, carverId, carver);
        Registry.register(BuiltinRegistries.CONFIGURED_CARVER, carverId, configuredCarver);
        return carverId;
    }

}
