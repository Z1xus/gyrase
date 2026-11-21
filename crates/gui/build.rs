use slint_build::{CompilerConfiguration, EmbedResourcesKind};

fn main() -> Result<(), slint_build::CompileError> {
	let config = CompilerConfiguration::new().embed_resources(EmbedResourcesKind::EmbedFiles);
	slint_build::compile_with_config("ui/app.slint", config)
}
