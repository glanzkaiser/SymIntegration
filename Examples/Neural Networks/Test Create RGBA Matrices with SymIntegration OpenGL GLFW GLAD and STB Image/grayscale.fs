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

	/* // Grayscale
	FragColor = texture(texture2, TexCoord);
	float average = (FragColor.r + FragColor.g + FragColor.b) / 3.0 ;
	FragColor = vec4(average, average, average, 1.0);
	*/

	
	// Weighted grayscale
	FragColor = texture(texture2, TexCoord);
	float average = 0.2126*FragColor.r + 0.7152*FragColor.g + 0.0722*FragColor.b ;
	FragColor = vec4(average, average, average, 1.0);
	 
}