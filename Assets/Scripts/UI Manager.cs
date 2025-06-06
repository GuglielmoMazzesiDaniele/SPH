using UnityEngine;
using UnityEngine.UIElements;

public class UIManager : MonoBehaviour
{
    public float smoothingFactor = 0.1f;
    
    private float _smoothedFPS;
    private Label _fpsLabel;

    private void Start()
    {
        // Extracting the UI Document
        var UI = GetComponent<UIDocument>();
        
        // Extracting the FPS Label
        _fpsLabel = UI.rootVisualElement.Q<Label>("FPS_Label");
    }

    void Update()
    {
        // Computing the current FPS
        var currentFrameFPS = 1.0f / Time.deltaTime;
        
        // Computing the smoothed FPS
        _smoothedFPS = Mathf.Lerp(_smoothedFPS, currentFrameFPS, smoothingFactor);
        
        // Updating the UI
        _fpsLabel.text = $"FPS: {_smoothedFPS:F1}";
    }
}
