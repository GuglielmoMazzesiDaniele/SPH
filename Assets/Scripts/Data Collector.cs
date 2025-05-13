using UnityEngine;

public class DataCollector : MonoBehaviour
{
    public float collectionDuration = 10f;   

    private bool _collecting = false;
    private float _timeElapsed = 0f;
    private float _totalFPS = 0f;
    private int _frameCount = 0;

    void Update()
    {
        // Start collecting
        if (!_collecting && Input.GetKeyDown(KeyCode.C))
        {
            _collecting = true;
            _timeElapsed = 0f;
            _totalFPS = 0f;
            _frameCount = 0;
            Debug.Log("Started FPS collection...");
        }

        // If not collecting, return
        if (!_collecting)
            return;
        
        // Computing current frame stats and adding to the pile
        var currentFPS = 1f / Time.deltaTime;
        _totalFPS += currentFPS;
        _frameCount++;
        _timeElapsed += Time.deltaTime;

        // If not enough time has passed, return
        if (_timeElapsed <= collectionDuration)
            return;
        
        // Ending data collection and printing
        var averageFPS = _totalFPS / _frameCount;
        Debug.Log($"[TimedFPSCollector] Average FPS over {collectionDuration} seconds: {averageFPS:F2}");
        _collecting = false;
    }
}
