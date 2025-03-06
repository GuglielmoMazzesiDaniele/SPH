using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.Rendering.Universal;
using UnityEngine.Rendering.RenderGraphModule;
using UnityEngine.Rendering.RenderGraphModule.Util;

// Heavily inspired by video "Introduction to the Render Graph in Unity 6" by Unity
// (https://www.youtube.com/watch?v=U8PygjYAF7A&ab_channel=Unity)

public class CustomVignetteRendererFeature : ScriptableRendererFeature
{
    public RenderPassEvent injectionPoint = RenderPassEvent.AfterRenderingPostProcessing;
    public Material material;
    [Range(0.0f, 1.0f)] public float intensity;
    [Range(0.0f, 1.0f)] public float smoothness;  
    
    private CustomVignetteRenderPass _customVignetteRenderPass;

    /// <inheritdoc/>
    public override void Create()
    {
        // Creating the scriptable pass injected in the rendering pipeline
        _customVignetteRenderPass = new CustomVignetteRenderPass
        {
            // Configures where the render pass should be injected.
            renderPassEvent = RenderPassEvent.AfterRenderingPostProcessing
        };
    }

    // Here you can inject one or multiple render passes in the renderer.
    // This method is called when setting up the renderer once per-camera.
    public override void AddRenderPasses(ScriptableRenderer renderer, ref RenderingData renderingData)
    {
        // If the material is not set, return
        if (material == null)
        {
            Debug.LogError("Material must be set in the Material Based Post Processing Feature.");
            return;
        }
        
        // Setting the pass material
        _customVignetteRenderPass.Setup(material, intensity, smoothness);
        
        // Injecting the pass in the rendering pipeline
        renderer.EnqueuePass(_customVignetteRenderPass);
    }
    
    private class CustomVignetteRenderPass : ScriptableRenderPass
    {
        // Pass related fields
        private const string PassName = "Custom Vignette";
        private Material _material;
        
        // Shader relates fields
        private readonly int _intensityID = Shader.PropertyToID("_Intensity");
        private readonly int _smoothnessID = Shader.PropertyToID("_Smoothness");
        
        public void Setup(Material material, float intensity, float smoothness)
        {
            // Setting the material
            _material = material;
            
            // Setting the material's variables
            _material.SetFloat(_intensityID, intensity);
            _material.SetFloat(_smoothnessID, smoothness);

            // The pass requires a temporary texture
            requiresIntermediateTexture = true;
        }

        // RecordRenderGraph is where the RenderGraph handle can be accessed, through which render passes can be added to the graph.
        // FrameData is a context container through which URP resources can be accessed and managed.
        public override void RecordRenderGraph(RenderGraph renderGraph, ContextContainer frameData)
        {
            // Extracting the resource data from the current stage of the pipeline
            var resourceData = frameData.Get<UniversalResourceData>();
            
            // Verifying if that current pass is writing directly to the screen back buffer
            if(resourceData.isActiveTargetBackBuffer)
                Debug.LogError($"Skipping render pass {PassName}, it requires an intermediate ColorTexture");
            
            // Extracting the current color texture
            var source = resourceData.activeColorTexture;
            
            // Extracting the texture description
            var destinationDescription = renderGraph.GetTextureDesc(source);
            
            // Fine-tuning the description
            destinationDescription.name = $"CameraColor-{PassName}";
            destinationDescription.clearBuffer = false;
            
            // Creating the destination texture
            var destination = renderGraph.CreateTexture(destinationDescription);
            
            // Initializing the parameters for the blit
            RenderGraphUtils.BlitMaterialParameters blitMaterialParameters = new (source, destination, _material, 0);
            
            // Executing the material on the destination texture
            renderGraph.AddBlitPass(blitMaterialParameters, PassName);
            
            // Overriding the current camera color with the just computed texture
            resourceData.cameraColor = destination;
        }
    }
}
