#version 330 core
out vec4 FragColor;

in vec3 ourColor;
in vec2 TexCoord;

// texture samplers
uniform sampler2D texture1;
uniform sampler2D texture2;
const float offset = 1.0 / 300.0;

void main()
{
	// linearly interpolate between both textures (10% texture1, 90% texture2)
	FragColor = mix(texture(texture1, TexCoord), texture(texture2, TexCoord), 0.9);

	// convolution kernel Sobel x
	vec2 offsets[9]= vec2[] (
	vec2(-offset, offset), // top-left
	vec2(0.0f, offset), // top-center
	vec2(offset, offset), // top-right
	vec2(-offset, 0.0f), // center-left
	vec2(0.0f, 0.0f), // center
	vec2(offset, 0.0f), //center-right
	vec2(-offset, -offset), // bottom-left
	vec2(0.0f, -offset), // bottom-center
	vec2(offset, -offset) // bottom-right
	);
	float kernel[9] = float[](
		-1, 0, 1,
		-2, 0, 2,
		-1, 0, 1	
	);
	vec3 sampleTex[9];
	for(int i = 0; i < 9; ++i)
	{
		sampleTex[i]= vec3(texture(texture2, TexCoord.st + offsets[i]));	
	}
	vec3 col = vec3(0.0);
	for(int i = 0; i < 9; ++i)
	{
		col += sampleTex[i] * kernel[i];
	}
	FragColor = vec4(col, 1.0);
	
}