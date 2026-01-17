"use client";

import { Canvas, useFrame } from "@react-three/fiber";
import { Stars, Sparkles, Environment, Float, Edges, Line } from "@react-three/drei";
import { Bloom, EffectComposer } from "@react-three/postprocessing";
import { useRef, useMemo, Suspense } from "react";
import * as THREE from "three";
import Link from "next/link";

// Custom Shader for the Glowy Rainbow Spectrum Cube
const RainbowSpectrumShader = {
  uniforms: {
    uTime: { value: 0 },
    uColor: { value: new THREE.Color("#050510") },
  },
  vertexShader: `
    varying vec2 vUv;
    varying vec3 vNormal;
    varying vec3 vPosition;
    void main() {
      vUv = uv;
      vNormal = normalize(normalMatrix * normal);
      vPosition = position;
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    }
  `,
  fragmentShader: `
    uniform float uTime;
    uniform vec3 uColor;
    varying vec2 vUv;
    varying vec3 vNormal;
    varying vec3 vPosition;

    // Simplex 2D noise
    vec3 permute(vec3 x) { return mod(((x*34.0)+1.0)*x, 289.0); }
    float snoise(vec2 v) {
      const vec4 C = vec4(0.211324865405187, 0.366025403784439, -0.577350269189626, 0.024390243902439);
      vec2 i  = floor(v + dot(v, C.yy) );
      vec2 x0 = v -   i + dot(i, C.xx);
      vec2 i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
      vec4 x12 = x0.xyxy + C.xxzz;
      x12.xy -= i1;
      i = mod(i, 289.0);
      vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 )) + i.x + vec3(0.0, i1.x, 1.0 ));
      vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy), dot(x12.zw,x12.zw)), 0.0);
      m = m*m ; m = m*m ;
      vec3 x = 2.0 * fract(p * C.www) - 1.0;
      vec3 h = abs(x) - 0.5;
      vec3 ox = floor(x + 0.5);
      vec3 a0 = x - ox;
      m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );
      vec3 g;
      g.x  = a0.x  * x0.x  + h.x  * x0.y;
      g.yz = a0.yz * x12.xz + h.yz * x12.yw;
      return 130.0 * dot(m, g);
    }

    vec3 hsv2rgb(vec3 c) {
      vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
      vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
      return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
    }

    void main() {
      vec2 uv = vUv;
      
      // Slow, heavy fluid noise
      float n1 = snoise(uv * 1.5 + uTime * 0.05);
      float n2 = snoise(uv * 3.0 - uTime * 0.03);
      
      // ABYSSAL SPECTRUM GHOST
      // Faint, spectral clouds that drift in the void
      float drift = n1 * n2;
      
      // Extremely low saturation (0.3) and very low brightness (0.15)
      // This makes the spectrum look like a ghost or a faint nebula reflection
      float hue = fract(uv.x * 0.2 + uv.y * 0.2 + n1 * 0.1 + uTime * 0.02);
      vec3 spectralFog = hsv2rgb(vec3(hue, 0.4, 0.12));
      
      // Near-zero base color
      vec3 darkVoid = vec3(0.001, 0.001, 0.002);
      
      // Fresnel for the edge definition (also very dim)
      float fresnel = pow(1.0 - abs(dot(vNormal, vec3(0,0,1))), 4.0);
      
      // Combine. The spectrum is barely there, appearing as "fog" in the blackness.
      vec3 finalColor = mix(darkVoid, spectralFog, smoothstep(-0.5, 1.0, drift) * 0.8 + fresnel * 0.2);
      
      // Subtly add back a tiny bit of rainbow to the Fresnel edges only
      finalColor += fresnel * spectralFog * 2.0;
      
      // Transparent but very dark
      gl_FragColor = vec4(finalColor, 0.7);
    }
  `,
  transparent: true,
  side: THREE.DoubleSide,
};

function FancyCube() {
  const meshRef = useRef<THREE.Mesh>(null);
  const material = useMemo(() => new THREE.ShaderMaterial(RainbowSpectrumShader), []);

  useFrame((state) => {
    if (meshRef.current) {
      meshRef.current.rotation.y = state.clock.elapsedTime * 0.25;
      meshRef.current.rotation.x = state.clock.elapsedTime * 0.1;
    }
    material.uniforms.uTime.value = state.clock.elapsedTime;
  });

  const size = 2.2;

  return (
    <group position={[2.5, 0, 0]}>
      <Float speed={2} rotationIntensity={0.5} floatIntensity={1}>
        <mesh ref={meshRef} material={material}>
          <boxGeometry args={[size, size, size]} />
          {/* GLOWY WHITE OUTLINES */}
          <Edges 
            threshold={15} 
            color="#ffffff" 
            scale={1.01}
          >
             <meshBasicMaterial color="#ffffff" toneMapped={false} />
          </Edges>
          
          {/* Inner structure for extra depth */}
          <Edges 
            threshold={15} 
            color="#ffffff" 
            scale={0.99}
          >
             <meshBasicMaterial color="#ffffff" transparent opacity={0.2} toneMapped={false} />
          </Edges>
        </mesh>
      </Float>
    </group>
  );
}

export default function Home() {
  return (
    <div className="flex min-h-screen w-full flex-col md:flex-row items-center justify-center p-8 md:p-24 relative overflow-hidden bg-black">
      
      <div className="absolute inset-0 z-0">
        <Canvas dpr={[1, 2]} camera={{ position: [0, 0, 8], fov: 40 }}>
          <ambientLight intensity={0.2} />
          <Stars radius={300} depth={60} count={20000} factor={8} saturation={0} fade speed={1} />
          <Sparkles scale={100} count={1000} size={2} speed={0.2} opacity={0.15} color="#ffffff" />
          
          <Suspense fallback={null}>
            <Environment resolution={1024}>
              <group rotation={[-Math.PI / 4, 0, 0]}>
                <mesh scale={200}>
                  <sphereGeometry args={[1, 64, 64]} />
                  <meshBasicMaterial color="#000000" side={THREE.BackSide} />
                </mesh>
              </group>
            </Environment>

            <FancyCube />

            <EffectComposer enableNormalPass={false}>
              <Bloom 
                luminanceThreshold={1} 
                mipmapBlur 
                intensity={1.2} // Increased intensity for glowy effect
                radius={0.4} 
              />
            </EffectComposer>
          </Suspense>
        </Canvas>
      </div>

      <div className="flex-1 z-10 flex flex-col items-center md:items-start text-center md:text-left space-y-6">
        <h1 className="text-6xl md:text-8xl font-bold tracking-tighter glow-text text-white">
          PROTEUS
        </h1>
        <p className="text-xl md:text-2xl text-gray-300 max-w-md font-light">
          Automated Molecular Dynamics Pipeline for Polymer Engineering.
        </p>
        <div className="flex gap-4 pt-4">
          <Link href="/simulation" className="px-8 py-3 bg-white text-black font-bold rounded-full hover:bg-gray-200 transition-all glow-box">
            Start Simulation
          </Link>
          <button className="px-8 py-3 border border-white text-white font-bold rounded-full hover:bg-white/10 transition-all">
            Documentation
          </button>
        </div>
      </div>

      <div className="flex-1 pointer-events-none" />
    </div>
  );
}
