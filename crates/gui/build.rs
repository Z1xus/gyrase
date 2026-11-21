use slint_build::{CompilerConfiguration, EmbedResourcesKind};

fn main() -> Result<(), Box<dyn std::error::Error>> {
	let config = CompilerConfiguration::new().embed_resources(EmbedResourcesKind::EmbedFiles);
	slint_build::compile_with_config("ui/app.slint", config)?;
	#[cfg(windows)]
	winresource::WindowsResource::new()
		.set_icon("assets/gyrase.ico")
		.compile()?;
	Ok(())
}
