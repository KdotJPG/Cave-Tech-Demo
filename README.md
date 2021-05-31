# K.jpg's Cave Tech Demo

A testbed for my various cave generation ideas. This is not intended to be used on playable Minecraft worlds in its current form.

Configure which generator is used, in [CaveTechDemo.java](https://github.com/KdotJPG/Cave-Tech-Demo/blob/master/src/main/java/jpg/k/cavetechdemo/CaveTechDemo.java).
Then run `gradlew runClient` in the main directory.

### Derivative Tunnel Closing Cave Generator

Computes the sum of squares of two seeds of OpenSimplex2 noise to produce smooth and roughly round tunnels.
Uses the cave's derivative vectors to check when caves are becoming less round and too distorted, to close off those
areas and turn them into dead ends. This addresses two potential stylistic blockers of noise tunnels at the same time.
Uses a cached column noise generator and a form of conditional noise layer skipping to enable fast generation.

![Derivative Tunnel Closing Caves 1](https://user-images.githubusercontent.com/8829856/120146071-83880100-c1b2-11eb-86b1-c933eab4b3be.png)

![Derivative Tunnel Closing Caves 2](https://user-images.githubusercontent.com/8829856/120146075-84209780-c1b2-11eb-9c05-8d422b74ad72.png)

[DerivativeTunnelClosingCaveCarver.java](https://github.com/KdotJPG/Cave-Tech-Demo/blob/master/src/main/java/jpg/k/cavetechdemo/carver/DerivativeTunnelClosingCaveCarver.java)

### Best Two Of Three Noises Cave Generator

This was a bit of a failed experiment. My goal was to smoothly and dynamically pick the best two of three noises
to carve caves with, but it didn't end up producing caves the way I intended. It's possible I made a math error,
but it's also possible the idea just doesn't work out as I had hoped. Leaving this in here in case anyone finds
parts of it useful. Particularly, the smooth weight selection formula might be great for other effects.
You can try it for yourself, but the results are a bit chaotic.

[BestTwoOfThreeNoisesCaveCarver.java](https://github.com/KdotJPG/Cave-Tech-Demo/blob/master/src/main/java/jpg/k/cavetechdemo/carver/BestTwoOfThreeNoisesCaveCarver.java)